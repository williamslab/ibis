#include <stdio.h>
#include <stdlib.h>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <genetio/marker.h>
#include <genetio/personbulk.h>
#include <genetio/personloopdata.h>
#include <genetio/personio.h>
#include <genetio/util.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <bitset>
#include <cstdlib>
#include <string.h>
#include <string>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <omp.h>
#include <sys/time.h>
#include <time.h>
#include <zlib.h>
#include <stdarg.h>
#include <map>




template<typename IO_TYPE>
class FileOrGZ {
	public:
		bool open(const char *filename, const char *mode);
		int getline();
		int printf(const char *format, ...);
		int close();

		static const int INIT_SIZE = 1024 * 50;

		// IO_TYPE is either FILE* or gzFile;
		IO_TYPE fp;

		// for I/O:
		char *buf;
		size_t buf_size;
		size_t buf_len;

	private:
		void alloc_buf();
};


template<typename IO_TYPE>
void FileOrGZ<IO_TYPE>::alloc_buf() {
	buf = (char *) malloc(INIT_SIZE);
	if (buf == NULL) {
		fprintf(stderr, "ERROR: out of memory!\n");
		exit(1);
	}
	buf_size = INIT_SIZE;
	buf_len = 0;
}

// open <filename> using standard FILE *
template<>
bool FileOrGZ<FILE *>::open(const char *filename, const char *mode) {
	// First allocate a buffer for I/O:
	alloc_buf();

	fp = fopen(filename, mode);
	if (!fp)
		return false;
	else
		return true;
}

// open <filename> as a gzipped file
template<>
bool FileOrGZ<gzFile>::open(const char *filename, const char *mode) {
	// First allocate a buffer for I/O:
	alloc_buf();

	fp = gzopen(filename, mode);
	if (!fp)
		return false;
	else
		return true;
}

template<>
int FileOrGZ<FILE *>::getline() {
	return ::getline(&buf, &buf_size, fp);
}

template<>
int FileOrGZ<gzFile>::getline() {
	int n_read = 0;
	int c;

	while ((c = gzgetc(fp)) != EOF) {
		// About to have read one more, so n_read + 1 needs to be less than *n.
		// Note that we use >= not > since we need one more space for '\0'
		if (n_read + 1 >= (int) buf_size) {
			const size_t GROW = 1024;
			char *tmp_buf = (char *) realloc(buf, buf_size + GROW);
			if (tmp_buf == NULL) {
				fprintf(stderr, "ERROR: out of memory!\n");
				exit(1);
			}
			buf_size += GROW;
			buf = tmp_buf;
		}
		buf[n_read] = (char) c;
		n_read++;
		if (c == '\n')
			break;
	}

	if (c == EOF && n_read == 0)
		return -1;

	buf[n_read] = '\0';

	return n_read;
}

template<>
int FileOrGZ<FILE *>::printf(const char *format, ...) {
	va_list args;
	int ret;
	va_start(args, format);

	ret = vfprintf(fp, format, args);
	va_end(args);
	return ret;
}

template<>
int FileOrGZ<gzFile>::printf(const char *format, ...) {
	va_list args;
	int ret;
	va_start(args, format);

	// NOTE: one can get automatic parallelization (in a second thread) for
	// gzipped output by opening a pipe to gzip (or bgzip). For example:
	//FILE *pipe = popen("gzip > output.vcf.gz", "w");
	// Can then fprintf(pipe, ...) as if it were a normal file.

	// gzvprintf() is slower than the code below that buffers the output.
	// Saw 13.4% speedup for processing a truncated VCF with ~50k lines and
	// 8955 samples.
	//  ret = gzvprintf(fp, format, args);
	ret = vsnprintf(buf + buf_len, buf_size - buf_len, format, args);
	if (ret < 0) {
		printf("ERROR: could not print\n");
		perror("printf");
		exit(10);
	}

	if (buf_len + ret > buf_size - 1) {
		// didn't fit the text in buf
		// first print what was in buf before the vsnprintf() call:
		gzwrite(fp, buf, buf_len);
		buf_len = 0;
		// now ensure that redoing vsnprintf() will fit in buf:
		if ((size_t) ret > buf_size - 1) {
			do { // find the buffer size that fits the last vsnprintf() call
				buf_size += INIT_SIZE;
			} while ((size_t) ret > buf_size - 1);
			free(buf);
			buf = (char *) malloc(buf_size);
		}
		// redo:
		ret = vsnprintf(buf + buf_len, buf_size - buf_len, format, args);
	}

	buf_len += ret;
	if (buf_len >= buf_size - 1024) { // within a tolerance of MAX_BUF?
		// flush:
		gzwrite(fp, buf, buf_len);
		buf_len = 0;
	}

	va_end(args);
	return ret;
}

template<>
int FileOrGZ<FILE *>::close() {
	assert(buf_len == 0);
	// should free buf, but I know the program is about to end, so won't
	return fclose(fp);
}

template<>
int FileOrGZ<gzFile>::close() {
	if (buf_len > 0)
		gzwrite(fp, buf, buf_len);
	return gzclose(fp);
}

//Class for storing the 3 types of bit sets required for window comparison.
class HomozygAndMiss{
	public:
		uint64_t hom1;//Set for homozygosity for 1st allele
		uint64_t hom2;//Set for homozygosity for 2nd allele.
		uint64_t miss;//Set for missing data
		HomozygAndMiss(){
			hom1 = 0;
			hom2 = 0;
			miss = 0;
		}
		void zero();//Currently unused?
		void addBits(uint64_t h1, uint64_t h2, uint64_t m);//setter function.
};



void HomozygAndMiss::zero(){
	hom1 = 0;
	hom2 = 0;
	miss = 0;
}

void HomozygAndMiss::addBits(uint64_t h1, uint64_t h2, uint64_t m){//Places single-bit values in the input into the bit sets.
	hom1+=h1;
	hom2+=h2;
	miss+=m;
}

HomozygAndMiss *seek(HomozygAndMiss *data, uint64_t indiv, uint64_t newMarkerBlocks){//Finds the starting location for a given individual in the data. Later HomozygAndMiss classes for that individual can be reached by adding 1 to the pointer.
	return data+(indiv*newMarkerBlocks);
}



//Representation of segment data for storage and processing.
class SegmentData{

	public:

		//NOTE: Stored segment endpoints correspond to the last error free window encountered.
		//Indexes of start and end markers.
		int startPos;
		int endPos;
		//Genetic positions of segment start and end.
		float startFloat;
		float endFloat;
		int errorCount; //Total errors considered part of the real ongoing segment.
		int trackedError; //Errors beyond the current end of the valid segment, extending into the ongoing consecutive error area.
		int trackedErrorIBD1; //For use in IBD2 segments to correctly calculate the errors in the following IBD1 segment when terminated.
		int ibdType; // Either 1 or 2 depending on nature of IBD for the segment.
		bool realIBD2; // Marks IBD1 segments that were contiguous with a stored real IBD2 segment.
		const char* chrom;
		SegmentData(){
			startPos = -1;
			startFloat = -1;//genetic start position of ongoing segment.
			errorCount = 0;//error count of ongoing segment.
			trackedError = 0;
			trackedErrorIBD1 = 0;
			endFloat = -1;
			endPos = -100;//position of last error free window. Used to backtrack segment endings to the last error free window.

			realIBD2 = false;
		}
		/*SegmentData(const SegmentData &oldSeg){
		  ind1 = oldSeg.ind1;
		  ind2 = oldSeg.ind2;
		  ibdType = oldSeg.ibdType;
		  startPos = -1;
		  startFloat = -1;//genetic start position of ongoing segment.
		  errorCount = 0;//error count of ongoing segment.
		  trackedError = 0;
		  trackedErrorIBD1 = 0;
		  endFloat = -1;
		  endPos = -100;//position of last error free window. Used to backtrack segment endings to the last error free window.
		  chrom = oldSeg.chrom;
		  realIBD2 = false;
		  }*/

		template<class IO_TYPE>	
			void printSegment(FileOrGZ<IO_TYPE> &pFile, int ind1, int ind2);

		float cMLength();
		int markerLength();
		bool errorCheck(int newError, int newEndPos, float errorRate);//Confirms segments error
		void updateSegmentEndpoints(int newEndPos, int newError, int newIBD1Error);
		void updateSegmentStartpoints(int newStartPos);
		bool checkSegment(float min_length, int marker_length);
		bool hasStart();
		void selfErase();
		void trackIBD1Errors();	


};


inline float SegmentData::cMLength(){return endFloat-startFloat;}
inline int SegmentData::markerLength(){return endPos-startPos;}


template<class IO_TYPE>
void SegmentData::printSegment(FileOrGZ<IO_TYPE> &pFile, int ind1, int ind2){

	pFile.printf("%s\t%s\t%s\t%i\t%i\tIBD%i\t%f\t%f\t%f\t%i\t%i\t%f\n", PersonLoopData::_allIndivs[ind1]->getId(),PersonLoopData::_allIndivs[ind2]->getId(),chrom, Marker::getMarker(startPos)->getPhysPos(), Marker::getMarker(endPos)->getPhysPos(),ibdType, startFloat, endFloat, cMLength(), markerLength(), int(errorCount), float(errorCount) / float(endPos - startPos));

}

//Return true if the segment with the new modifications fails the error threshold test: true if fail, false if pass.
inline bool SegmentData::errorCheck(int newError, int newEndPos, float errorRate){
	return ((errorCount + trackedError + newError) / (newEndPos - startPos)) > errorRate;
}

void SegmentData::updateSegmentEndpoints(int newEndPos, int newError, int newIBD1Error){
	if(newError<=0){
		endPos = newEndPos;
		endFloat = Marker::getMarker(endPos)->getMapPos();
		errorCount+=trackedError;
		trackedError = 0;
		trackedErrorIBD1 = 0;
	}
	else{//Won't update real endpoint if there was an error in the window, but tracks the information anyway for error analyses.
		trackedError += newError;
		trackedErrorIBD1 += newIBD1Error;
	}
}

void SegmentData::updateSegmentStartpoints(int newStartPos){
	startPos = newStartPos;
	startFloat = Marker::getMarker(startPos)->getMapPos();
	trackedError = 0;
	trackedErrorIBD1 = 0;
	errorCount = 0;
	realIBD2 = false;
}
//Confirms if segment meets relevant conditions to be a valid segment.
inline bool SegmentData::checkSegment(float min_length, int marker_length){
	return ((endPos > startPos) && (startPos >= 0) && (realIBD2 || (((endPos - startPos) > marker_length) && ((endFloat - startFloat) > min_length))));
}
//Invalid segments are reprsented by a startPos of -1.
inline bool SegmentData::hasStart(){
	return startPos >= 0;
}

void SegmentData::selfErase(){//TODO This may be streamlinable. Will check later.
	startPos = -1;
	startFloat = -1;
	endPos = -100;
	endFloat = -1;
	errorCount = 0;
	realIBD2 = false;
}


void storeSegment(SegmentData& currentSegment, std::vector<SegmentData> &storedSegments, float &totalIBD1, float &totalIBD2){
	storedSegments.push_back((const SegmentData)(currentSegment));
	totalIBD1+=currentSegment.cMLength() * (currentSegment.ibdType == 1);
	totalIBD2+=currentSegment.cMLength() * (currentSegment.ibdType == 2);
	currentSegment.selfErase();
}

//Find the lowest order non-zero bit in the uint64_t, and return it in a uint64_t with only that bit set to 1.
inline uint64_t lowestSetBit(uint64_t bitSet)
{
	return (-bitSet) & bitSet;
}

//Determine if the two given segments overlap
inline bool segmentOverlap(int start1, int start2, int end1, int end2){
	return (end2 >= start1) && (end1 >= start2);//(start1 <= end1 && start1 >= start2) || (end1<=end2 && end1>=start2);

}



//Segment termination code.
//Tries to merge given segment with the stored "last" segment. If that fails, checks the ongoing segment against all the criteria for segment validity. If it succeeds, it prints the stored segment and stores the ongoing one in the "last" segment slot.
//Returns a boolean describing if the segment in storage is a valid segment according to the thresholds, which, as a side effect, will also always be true if the segment just analyzed was valid.
//template<class IO_TYPE>
bool endSegment(SegmentData& currentSegment, float min_length, int min_markers,  std::vector<SegmentData> &storedSegments, float &totalIBD1, float &totalIBD2){ 
	bool passedChecks = false;
	if(currentSegment.checkSegment(min_length, min_markers)){
		//currentSegment.printSegment(pFile);
		storeSegment(currentSegment, storedSegments, totalIBD1, totalIBD2);
		passedChecks = true;
	}
	//currentSegment.selfErase();
	return passedChecks;


}

//Breaks and prints underlying IBD1 segment after IBD2 segment confirmed to be real. Also updates IBD1 ongoing segment data.
//template<class IO_TYPE>
//template<class IO_TYPE>
bool handleIBD1PostIBD2(SegmentData &ibd1Segment,  std::vector<SegmentData> &storedSegments, float &totalIBD1, float &totalIBD2){
	SegmentData ibd2Segment = storedSegments.back();
	if(ibd2Segment.startPos > ibd1Segment.startPos && ibd1Segment.startPos >= 0){
		ibd1Segment.updateSegmentEndpoints(ibd2Segment.startPos,0,0);
		ibd1Segment.errorCount = -1;
		int ibd1EndPosHolder = ibd1Segment.endPos;
		//ibd1Segment.printSegment(pFile);
		storeSegment(ibd1Segment, storedSegments, totalIBD1, totalIBD2);
		ibd1Segment.updateSegmentStartpoints(ibd2Segment.endPos);
		ibd1Segment.updateSegmentEndpoints(ibd1EndPosHolder, 0, 0);
		ibd1Segment.errorCount = ibd2Segment.trackedErrorIBD1;//Necessary to put real error count back into ongoing IBD1 segment.
		ibd1Segment.realIBD2 = true;
		return true;
	}
	ibd1Segment.updateSegmentStartpoints(ibd2Segment.endPos);
	ibd1Segment.realIBD2 = true;
	return false;
}





//Determine the corresponding block of 64 and bit within that block for genotype data.
void findIndexAndShift(uint64_t marker, uint64_t indiv, uint64_t markerWindows, uint64_t prevChromBlocks, uint64_t &newIndex, uint64_t &newShift){
	newIndex=indiv * (markerWindows)+prevChromBlocks+((marker / 64));//Finds the 64-bit region in the transposed data where the value (first value) fits. Second goes in the next block
	newShift = marker % 64;//Finds how deep into that 64-bit regions the new bits need to be placed.
}

//Determine the corresponding block of 64 and bit within that block for missing genotype data in a difference sized matrix.
void findIndexAndShiftMiss(uint64_t marker, uint64_t indiv, uint64_t newMarkerBlocks, uint64_t &newIndex, uint64_t &newShift){
	newIndex=indiv * (newMarkerBlocks / 2) + (marker / 64);
	newShift = marker % 64;
}

//Find block corresponding to an individual and set of markers.
inline uint64_t findIndex(uint64_t markerWindow, uint64_t indiv, uint64_t newMarkerBlocks){
	return indiv*(2*newMarkerBlocks)+(2*markerWindow);
}

//Find blockcorresponding to an individual and set of markers in the missing data.
inline uint64_t findMissIndex(uint64_t markerWindow, uint64_t indiv, uint64_t newMarkerBlocks){
	return indiv*newMarkerBlocks+markerWindow;
}



//populate the IBD boolean and error boolean int with true or false depending on if there is valid IBD1 in the current window and if it had any errors in it.
//Populate the errorCount value with the number of errors. Necessary for merge, error density analyses, and skipping IBD2 analysis when unneeded. 
void compareWindowsIBD1(uint64_t hom11, uint64_t hom12, uint64_t hom21, uint64_t hom22, int &errorCount, bool &IBD){
	uint64_t fullCompare = 	(hom11 & hom22) | (hom12 & hom21);//The bitwise comparison for iBD1

	errorCount = __builtin_popcountl(fullCompare);
	IBD = errorCount < 2;
	return;
}

//populate the IBD boolean and error boolean int with true or false depending on if there is valid IBD1 in the current window and if it had any errors in it.
//Populate the errorCount value with the number of errors. Necessary for merge and error density analyses.
void compareWindowsIBD2(uint64_t hom11, uint64_t hom12, uint64_t hom21, uint64_t hom22, uint64_t missing1, uint64_t missing2, int &errorCount, bool &IBD){
	uint64_t fullCompare = (~(missing1 | missing2)) & ((hom11 ^ hom21) | (hom12 ^ hom22));//The bitwise comparison for IBD2.


	errorCount = __builtin_popcountl(fullCompare);
	IBD = errorCount < 3;
	return;
}

//Transpose the input genotypes into a matrix of individuals as rows and blocks of markers as columns to speed up our analyses. Populates transposedData with homozygosity bit sets, and missingData with bit sets of where data is missing.
void interleaveAndConvertData(HomozygAndMiss transposedData[], int64_t indBlocks64, uint64_t markerWindows, uint64_t markers, int64_t individuals, std::vector<int> &starts, std::vector<int> &ends,  uint64_t chromWindowStarts[]) {

	printf("Organizing genotype data for analysis... "); 
	fflush(stdout);
	uint64_t dataBlocks = ((individuals+3)/4);//Breakdown of the 8 bit block counts based on the number of individuals in the input. Not the blocks in our own internal data format.
	//Need to check later if redundant with blocks. These were different quantities before.

	uint64_t mask1 = 1;


	//std::fill(transposedData, transposedData + individuals * newMarkerBlocks, 0 );
	//std::fill(missingData, missingData + individuals * (newMarkerBlocks / 2), 0);


	uint8_t *data = new uint8_t[dataBlocks];

	for (int chrIndex = 0; chrIndex < Marker::getNumChroms(); chrIndex++){
		for (int m = 0; m < Marker::getNumChromMarkers(chrIndex); m++){
			PersonIO<PersonLoopData>::readGenoRow(data, (int)(dataBlocks));//Reading in one row of the .bed file input to be transposed.
			//printf("\nchrom%s %i %u\n",Marker::getChromName(chrIndex), m, data[0]);
			uint64_t *data64 = (uint64_t *)data;//Reformat to allow blocks of 64 bits to be extracted and processed at once from an otherwise 8 bit format.

			uint8_t *data8 = data;
			for (int64_t b = 0; b < indBlocks64; b++) {

				uint64_t individuals1;
				uint64_t individuals2;
				if (b < (int64_t)(individuals) / 64) {
					individuals1 = data64[b * 2];//First set of 32 2-bit individuals to be combined into two blocks of 64 later.
					individuals2 = data64[b * 2 + 1];//Second set of 32 2-bit individuals to be combined into two blocks of 64 later. TODO May not be a speedup in new data ordering.
				}
				else {
					//Handle if the first set of individuals does not fill out a set of 64 completely and needs to be supplemented with 0s.
					if(individuals%64<32){
						individuals1 = 0;
						individuals2 = 0;
						int count = 0;
						for (int64_t x = ((individuals-1) / 4); x>=b*16; x--) {
							individuals1 = individuals1 << 8;
							individuals1 += data8[x];
							count++;
						}
					}
					else{//Handle if the second set of individuals does not fill out a set of 64 completely and needs to be supplemented with 0s.
						individuals1 = data64[b * 2];
						individuals2 = 0;
						int count = 0;
						for (int64_t x = ((individuals-1) / 4); x>=b*16+8; x--) {
							individuals2 = individuals2 << 8;
							individuals2 += data8[x];

							count++;

						}
					}

				}
				uint64_t newShift;
				uint64_t newIndex;

				//converts 2-bit format to 2 1-bit homozygotes and places single bit value into corresponding point in transposed data
				//Leaving this in two pieces for now, as there are differences due to how the two sets of individuals are placed in the final set of 64 index wise. Will combine later..
				for(uint64_t indivShift=0; indivShift<32; indivShift++){
					uint64_t indiv = indivShift+(b * 64);
					if((int64_t)(indiv)>=individuals)
						break;


					//Put the two relevant bits for the current genotype at the end of the set and mask off the rest.
					uint64_t bit1 = (individuals1 >> (indivShift * 2)) & mask1;
					uint64_t bit2 = ((individuals1 >> (indivShift * 2+1)) & mask1);


					//Create homozygous set data.
					uint64_t hom1 = bit1 & bit2;
					uint64_t hom2 = ((~bit1) & (~bit2)) & mask1;//TODO POSSIBLE SUBOPTIMAL
					//Create missing set data
					uint64_t missBits = bit1 & ~bit2 & mask1;

					//Find and emplace data in transposed matrix.
					findIndexAndShift(m, indiv, markerWindows, chromWindowStarts[chrIndex], newIndex, newShift);
					transposedData[newIndex].addBits(hom1<<newShift,hom2<<newShift,missBits<<newShift);

					//uint64_t hom1Shift = hom1<<newShift;
					//transposedData[newIndex]+=hom1Shift;

					//uint64_t hom2Shift = hom2<<newShift;
					//transposedData[newIndex+1]+=hom2Shift;



					//findIndexAndShiftMiss(m, indiv, newMarkerBlocks, newIndex, newShift);
					//uint64_t missShift = missBits<<newShift;
					//missingData[newIndex]+=missShift;

				}
				for(uint64_t indivShift=32; indivShift<64; indivShift++){
					uint64_t indiv = indivShift+(b * 64);
					if((int64_t)(indiv)>=individuals)
						break;
					uint64_t bit1 = (individuals2 >> (indivShift * 2)) & mask1;
					uint64_t bit2 = ((individuals2 >> (indivShift * 2+1)) & mask1);

					uint64_t hom1 = bit1 & bit2;
					uint64_t hom2 = ((~bit1) & (~bit2)) & mask1;

					uint64_t missBits = bit1 & ~bit2 & mask1;

					findIndexAndShift(m, indiv, markerWindows, chromWindowStarts[chrIndex], newIndex, newShift);
					transposedData[newIndex].addBits(hom1<<newShift,hom2<<newShift,missBits<<newShift);

					//uint64_t hom1Shift = hom1<<newShift;
					//transposedData[newIndex]+=hom1Shift;


					//uint64_t hom2Shift = hom2<<newShift;
					//transposedData[newIndex+1]+=hom2Shift;
					//findIndexAndShiftMiss(m, indiv, newMarkerBlocks, newIndex, newShift);
					//uint64_t missShift = missBits<<newShift;
					//missingData[newIndex]+=missShift;

				}

			}
		}
	}
	delete[] data;
	PersonIO<PersonLoopData>::closeGeno();
	printf("done.\n");
}

template<class IO_TYPE>
void internalIteration(HomozygAndMiss transposedData[], int numIndivs, int indBlocks, int numMarkers, int markerWindows, std::vector<int> &starts, std::vector<int> &ends,  std::string filename, std::string extension, uint64_t chromWindowBoundaries[][2], uint64_t chromWindowStarts[], int min_markers, float min_length, int min_markers2, float min_length2, float errorThreshold, float errorThreshold2, bool ibd2, bool printCoef, int numThreads, float min_coef, float fudgeFactor, int indiv1, int indiv2, FileOrGZ<IO_TYPE> pFileArray[], FileOrGZ<IO_TYPE> classFileArray[], std::vector<SegmentData> storedSegArray[], float ibd1Thresholds[], float ibd2Thresholds[], float totalGeneticLength){


	FileOrGZ<IO_TYPE> pFile = pFileArray[0];
	FileOrGZ<IO_TYPE> classFile = classFileArray[0];
	std::vector<SegmentData> storedSegs = storedSegArray[0];
	float totalIBD1Length = 0;
	float totalIBD2Length = 0;

	SegmentData currentIBD1Segment;// = new SegmentData;
	SegmentData currentIBD2Segment;// = new SegmentData;
	currentIBD1Segment.ibdType = 1;


	//IBD2 equivalents of the previously described variables.
	currentIBD2Segment.ibdType = 2;

	HomozygAndMiss *ind1Data = seek(transposedData,indiv1,markerWindows);//Stores the indices of the startpoints for the two individuals in the stored homozygous data and missing data.
	HomozygAndMiss *ind2Data = seek(transposedData,indiv2,markerWindows);


	for(int chr = 0; chr < Marker::getNumChroms(); chr++){
		currentIBD1Segment.chrom = Marker::getChromName(chr);
		currentIBD2Segment.chrom = Marker::getChromName(chr);
		for(uint64_t markerWindow=chromWindowBoundaries[chr][0]; markerWindow<=chromWindowBoundaries[chr][1]; markerWindow++){
			//Booleans for the status of the IBD conclusions for the current window comparison.
			bool isIBD1;
			bool isIBD2;

			int windowErrorCount;
			int windowErrorCount2;


			compareWindowsIBD1(ind1Data->hom1, ind1Data->hom2, ind2Data->hom1, ind2Data->hom2, windowErrorCount, isIBD1);

			//Don't bother checking IBD2 if there are too many errors in IBD1.
			if(ibd2 && isIBD1){
				compareWindowsIBD2(ind1Data->hom1, ind1Data->hom2, ind2Data->hom1, ind2Data->hom2, ind1Data->miss, ind2Data->miss, windowErrorCount2, isIBD2);
			}
			else{
				isIBD2=false;
				windowErrorCount2=3;
			}

			if(ibd2){//The check for if IBD2 is enabled
				if(isIBD2){
					int currentEndPos = ends[markerWindow];
					if(windowErrorCount2>0){
						if(currentIBD2Segment.hasStart() && currentIBD2Segment.errorCheck(windowErrorCount2, currentEndPos, errorThreshold2)){//Update potential endpoint if error threshold passed.
							bool realSeg = endSegment(currentIBD2Segment, min_length2, min_markers2,  storedSegs, totalIBD1Length, totalIBD2Length);
							if(realSeg){
								handleIBD1PostIBD2(currentIBD1Segment,  storedSegs, totalIBD1Length, totalIBD2Length);
							}
							currentIBD2Segment.selfErase();//May be unnecessary relic from earlier draft. Some of these are, some aren't.
						}
						else{//Update endpoints if no errors.
							currentIBD2Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount2, windowErrorCount);
						}

					}
					else{

						if(!currentIBD2Segment.hasStart())//Update startopint if there is no current real segment.
							currentIBD2Segment.updateSegmentStartpoints(starts[markerWindow]);
						currentIBD2Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount2, windowErrorCount);//Update endpoint always.
					}

				}
				else {//If no IBD, attempt to store segment if valid.
					bool realSeg = endSegment(currentIBD2Segment, min_length2, min_markers2,  storedSegs, totalIBD1Length, totalIBD2Length);
					if(realSeg){

						handleIBD1PostIBD2(currentIBD1Segment,  storedSegs, totalIBD1Length, totalIBD2Length);
					}
					currentIBD2Segment.selfErase();
				}
			}




			if(isIBD1){
				int currentEndPos = ends[markerWindow];
				if(windowErrorCount>0){
					if(currentIBD1Segment.hasStart() && currentIBD1Segment.errorCheck(windowErrorCount, currentEndPos, errorThreshold)){
						endSegment(currentIBD1Segment, min_length, min_markers,  storedSegs, totalIBD1Length, totalIBD2Length);
						currentIBD1Segment.selfErase();

					}
					else{
						currentIBD1Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount, windowErrorCount);
					}

				}
				else{
					if(!currentIBD1Segment.hasStart())
						currentIBD1Segment.updateSegmentStartpoints(starts[markerWindow]);
					currentIBD1Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount, windowErrorCount);
				}
			}
			else {
				endSegment(currentIBD1Segment, min_length, min_markers,  storedSegs, totalIBD1Length, totalIBD2Length);
				currentIBD1Segment.selfErase();
			}

			ind1Data++;
			ind2Data++;

		}
		//Handle the end of chromosomes, as there is no end of IBD or error density increase to trigger these segments to print otherwise.

		if(ibd2){
			//Attempt to print IBD2 at the end of the chromosome if there is a valid segment.
			bool realSeg = endSegment(currentIBD2Segment, min_length2, min_markers2,  storedSegs, totalIBD1Length, totalIBD2Length);
			if(realSeg){
				handleIBD1PostIBD2(currentIBD1Segment,  storedSegs, totalIBD1Length, totalIBD2Length);
			}
		}
		//Attempt to print IBD1 at the end of the chromosome if there is a valid segment.
		endSegment(currentIBD1Segment, min_length, min_markers,  storedSegs, totalIBD1Length, totalIBD2Length);
		currentIBD1Segment.selfErase();
		currentIBD2Segment.selfErase();
	}



	//Perform relationship inference and compare to coefficient threshold.
	float relCoef = (0.25*totalIBD1Length+0.5*totalIBD2Length)/totalGeneticLength + fudgeFactor;//NOTE: Currently uses IBD2-based coefficient threshold whether or not IBD2 is enabled.
	float ibd1TotalMod = (1 - (totalIBD1Length + 4 * fudgeFactor)/totalGeneticLength);

	if(relCoef >= min_coef){
		if(printCoef){
			int classVal;
			if(ibd2){
				for(classVal=0; classVal<8; classVal++){
					if((relCoef)>ibd2Thresholds[classVal]){
						break;
					}
				}
			}
			else{
				for(classVal=0; classVal<8; classVal++){
					if(ibd1TotalMod<ibd1Thresholds[classVal]){
						break;
					}
				}
			}
			if(classVal==8)
				classVal=-1;
			//omp_set_lock(&lock);	
			classFile.printf("%s %s %f %f %i %i\n",PersonLoopData::_allIndivs[indiv1]->getId(),PersonLoopData::_allIndivs[indiv2]->getId(), relCoef, storedSegs.size(), totalIBD2Length/totalGeneticLength, classVal);
		}
		//omp_unset_lock(&lock);
		for (SegmentData seg: storedSegs){
			seg.printSegment(pFile,indiv1,indiv2);
		}
	}
	storedSegs.clear();
	//storedSegs.resize(maxSize);
	//delete currentIBD1Segment;
	//delete currentIBD2Segment;



}
//Loop over the individuals and windows and perform the segment analysis.
template<class IO_TYPE>
void segmentAnalysisMonoThread(HomozygAndMiss transposedData[], int numIndivs, int indBlocks, int numMarkers, int markerWindows, std::vector<int> &starts, std::vector<int> &ends,  std::string filename, std::string extension, uint64_t chromWindowBoundaries[][2], uint64_t chromWindowStarts[], int min_markers, float min_length, int min_markers2, float min_length2, float errorThreshold, float errorThreshold2, bool ibd2, bool printCoef, int numThreads, float min_coef, float fudgeFactor){	

	//omp_lock_t lock;//Lock for the coefficient/relationship class file
	//omp_init_lock(&lock);



	float totalGeneticLength = 0.0;//Value to use as whole input genetic length for calculating coefficients.
	for(int chr = 0; chr < Marker::getNumChroms(); chr++)
		totalGeneticLength += Marker::getMarker( Marker::getLastMarkerNum( chr ) )->getMapPos();
	float ibd1Thresholds[9];//Thresholds for class calling
	float ibd2Thresholds[9];
	ibd1Thresholds[0]=0;
	ibd1Thresholds[1]=0.365;//Hardcode for IBD1 only 1st threshold due to not following the formula.
	for(int i = 2; i<=8; i++)
		ibd1Thresholds[i] = 1-(1.0/pow(2.0,((i-1)+1.5)));

	//	{0.365, 1-(1.0/(2.0^(0+1.5))),1-(1.0/(2.0^(1+1.5))),1-(1.0/(2.0^(2+1.5))),1-(1.0/(2.0^(3+1.5))),1-(1.0/(2.0^(4+1.5))),1-(1.0/(2.0^(5+1.5))),1-(1.0/(2.0^(6+1.5))),1-(1.0/(2.0^(7+1.5)))};

	ibd2Thresholds[0]=0.475;
	for(int i = 1; i<=8; i++)
		ibd2Thresholds[i] = (1.0/(pow(2.0,(i+1.5))));

	//Transposes the input matrix and handles the 8->64 bit format conversion..
	interleaveAndConvertData(transposedData, indBlocks, markerWindows, numMarkers, numIndivs, starts, ends,  chromWindowStarts);

	//omp_set_dynamic(0);     // Explicitly disable dynamic teams
	//omp_set_num_threads(1); //Forces input specified thread number with default 1. 

	printf("Beginning segment detection with %i thread(s)...",numThreads);
	fflush(stdout);

	/*FileOrGZ<IO_TYPE> pFileArray[numThreads];
	  FileOrGZ<IO_TYPE> classFileArray[numThreads];
	  std::vector<SegmentData> storedSegArray[numThreads];


	  for(int thread = 0; thread < numThreads; thread++){	
	  std::string threadname;
	  threadname = filename +".seg."+std::to_string((1+thread))+extension;
	  bool success = pFileArray[thread].open(threadname.c_str(), "w");
	  if(!success){
	  printf("\nERROR: could not open output VCF file %s!\n", threadname.c_str());
	  perror("open");
	  exit(1);
	  }


	  if(printCoef){
	  std::string classFilename;
	  classFilename = filename +".coef."+std::to_string((1+thread))+extension;
	  success = classFileArray[thread].open(classFilename.c_str(), "w");
	  if(!success){
	  printf("\nERROR: could not open output VCF file %s!\n", classFilename.c_str());
	  perror("open");
	  exit(1);
	  }


	  classFileArray[thread].printf("Individual1\tIndividual2\tKinship_Coefficient\tIBD2_Fraction\tSegment_Count\tDegree\n");
	  }
	  }*/



	std::string threadname;
	FileOrGZ<IO_TYPE> pFile;
	FileOrGZ<IO_TYPE> classFile;
	threadname = filename +".seg"+extension;
	bool success = pFile.open(threadname.c_str(), "w");
	if(!success){
		printf("\nERROR: could not open output VCF file %s!\n", threadname.c_str());
		perror("open");
		exit(1);
	}
	if(printCoef){
		std::string classFilename;
		classFilename = filename +".coef."+extension;
		success = classFile.open(classFilename.c_str(), "w");
		if(!success){
			printf("\nERROR: could not open output VCF file %s!\n", classFilename.c_str());
			perror("open");
			exit(1);
		}


		classFile.printf("Individual1\tIndividual2\tKinship_Coefficient\tIBD2_Fraction\tSegment_Count\tDegree\n");
	}

	std::vector<SegmentData> storedSegs;
	for(uint64_t indiv1=0; indiv1<numIndivs; indiv1++){
		for(uint64_t indiv2=indiv1+1; indiv2<numIndivs; indiv2++){
			//FileOrGZ<IO_TYPE> pFile = pFileArray[0];
			//FileOrGZ<IO_TYPE> classFile = classFileArray[0];
			//std::vector<SegmentData> storedSegs = storedSegArray[0];
			float totalIBD1Length = 0;
			float totalIBD2Length = 0;

			SegmentData currentIBD1Segment;// = new SegmentData;
			SegmentData currentIBD2Segment;// = new SegmentData;
			currentIBD1Segment.ibdType = 1;


			//IBD2 equivalents of the previously described variables.
			currentIBD2Segment.ibdType = 2;

			HomozygAndMiss *ind1Data = seek(transposedData,indiv1,markerWindows);//Stores the indices of the startpoints for the two individuals in the stored homozygous data and missing data.
			HomozygAndMiss *ind2Data = seek(transposedData,indiv2,markerWindows);


			for(int chr = 0; chr < Marker::getNumChroms(); chr++){
				currentIBD1Segment.chrom = Marker::getChromName(chr);
				currentIBD2Segment.chrom = Marker::getChromName(chr);
				for(uint64_t markerWindow=chromWindowBoundaries[chr][0]; markerWindow<=chromWindowBoundaries[chr][1]; markerWindow++){
					//Booleans for the status of the IBD conclusions for the current window comparison.
					bool isIBD1;
					bool isIBD2;

					int windowErrorCount;
					int windowErrorCount2;


					compareWindowsIBD1(ind1Data->hom1, ind1Data->hom2, ind2Data->hom1, ind2Data->hom2, windowErrorCount, isIBD1);

					//Don't bother checking IBD2 if there are too many errors in IBD1.
					if(ibd2 && isIBD1){
						compareWindowsIBD2(ind1Data->hom1, ind1Data->hom2, ind2Data->hom1, ind2Data->hom2, ind1Data->miss, ind2Data->miss, windowErrorCount2, isIBD2);
					}
					else{
						isIBD2=false;
						windowErrorCount2=3;
					}

					if(ibd2){//The check for if IBD2 is enabled
						if(isIBD2){
							int currentEndPos = ends[markerWindow];
							if(windowErrorCount2>0){
								if(currentIBD2Segment.hasStart() && currentIBD2Segment.errorCheck(windowErrorCount2, currentEndPos, errorThreshold2)){//Update potential endpoint if error threshold passed.
									bool realSeg = endSegment(currentIBD2Segment, min_length2, min_markers2,  storedSegs, totalIBD1Length, totalIBD2Length);
									if(realSeg){
										handleIBD1PostIBD2(currentIBD1Segment,  storedSegs, totalIBD1Length, totalIBD2Length);
									}
									currentIBD2Segment.selfErase();//May be unnecessary relic from earlier draft. Some of these are, some aren't.
								}
								else{//Update endpoints if no errors.
									currentIBD2Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount2, windowErrorCount);
								}

							}
							else{

								if(!currentIBD2Segment.hasStart())//Update startopint if there is no current real segment.
									currentIBD2Segment.updateSegmentStartpoints(starts[markerWindow]);
								currentIBD2Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount2, windowErrorCount);//Update endpoint always.
							}

						}
						else {//If no IBD, attempt to store segment if valid.
							bool realSeg = endSegment(currentIBD2Segment, min_length2, min_markers2,  storedSegs, totalIBD1Length, totalIBD2Length);
							if(realSeg){

								handleIBD1PostIBD2(currentIBD1Segment,  storedSegs, totalIBD1Length, totalIBD2Length);
							}
							currentIBD2Segment.selfErase();
						}
					}




					if(isIBD1){
						int currentEndPos = ends[markerWindow];
						if(windowErrorCount>0){
							if(currentIBD1Segment.hasStart() && currentIBD1Segment.errorCheck(windowErrorCount, currentEndPos, errorThreshold)){
								endSegment(currentIBD1Segment, min_length, min_markers,  storedSegs, totalIBD1Length, totalIBD2Length);
								currentIBD1Segment.selfErase();

							}
							else{
								currentIBD1Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount, windowErrorCount);
							}

						}
						else{
							if(!currentIBD1Segment.hasStart())
								currentIBD1Segment.updateSegmentStartpoints(starts[markerWindow]);
							currentIBD1Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount, windowErrorCount);
						}
					}
					else {
						endSegment(currentIBD1Segment, min_length, min_markers,  storedSegs, totalIBD1Length, totalIBD2Length);
						currentIBD1Segment.selfErase();
					}

					ind1Data++;
					ind2Data++;

				}
				//Handle the end of chromosomes, as there is no end of IBD or error density increase to trigger these segments to print otherwise.

				if(ibd2){
					//Attempt to print IBD2 at the end of the chromosome if there is a valid segment.
					bool realSeg = endSegment(currentIBD2Segment, min_length2, min_markers2,  storedSegs, totalIBD1Length, totalIBD2Length);
					if(realSeg){
						handleIBD1PostIBD2(currentIBD1Segment,  storedSegs, totalIBD1Length, totalIBD2Length);
					}
				}
				//Attempt to print IBD1 at the end of the chromosome if there is a valid segment.
				endSegment(currentIBD1Segment, min_length, min_markers,  storedSegs, totalIBD1Length, totalIBD2Length);
				currentIBD1Segment.selfErase();
				currentIBD2Segment.selfErase();
			}



			//Perform relationship inference and compare to coefficient threshold.
			float relCoef = (0.25*totalIBD1Length+0.5*totalIBD2Length)/totalGeneticLength + fudgeFactor;//NOTE: Currently uses IBD2-based coefficient threshold whether or not IBD2 is enabled.
			float ibd1TotalMod = (1 - (totalIBD1Length + 4 * fudgeFactor)/totalGeneticLength);

			if(relCoef >= min_coef){
				if(printCoef){
					int classVal;
					if(ibd2){
						for(classVal=0; classVal<8; classVal++){
							if((relCoef)>ibd2Thresholds[classVal]){
								break;
							}
						}
					}
					else{
						for(classVal=0; classVal<8; classVal++){
							if(ibd1TotalMod<ibd1Thresholds[classVal]){
								break;
							}
						}
					}
					if(classVal==8)
						classVal=-1;
					//omp_set_lock(&lock);	
					classFile.printf("%s %s %f %f %i %i\n",PersonLoopData::_allIndivs[indiv1]->getId(),PersonLoopData::_allIndivs[indiv2]->getId(), relCoef, storedSegs.size(), totalIBD2Length/totalGeneticLength, classVal);
				}
				//omp_unset_lock(&lock);
				for (SegmentData seg: storedSegs){
					seg.printSegment(pFile,indiv1,indiv2);
				}
			}
			storedSegs.clear();
			//storedSegs.resize(maxSize);
			//delete currentIBD1Segment;
			//delete currentIBD2Segment;

			//internalIteration(transposedData, numIndivs, indBlocks, numMarkers, markerWindows, starts, ends,  filename, extension, chromWindowBoundaries, chromWindowStarts, min_markers, min_length, min_markers2, min_length2, errorThreshold, errorThreshold2, ibd2, printCoef, numThreads, min_coef, fudgeFactor, indiv1, indiv2,  pFileArray, classFileArray, storedSegArray, ibd1Thresholds, ibd2Thresholds, totalGeneticLength);
		}
	}

	//omp_destroy_lock(&lock);
	printf("done.\n");
}


//Loop over the individuals and windows and perform the segment analysis.
template<class IO_TYPE>
void segmentAnalysis(HomozygAndMiss transposedData[], int numIndivs, int indBlocks, int numMarkers, int markerWindows, std::vector<int> &starts, std::vector<int> &ends,  std::string filename, std::string extension, uint64_t chromWindowBoundaries[][2], uint64_t chromWindowStarts[], int min_markers, float min_length, int min_markers2, float min_length2, float errorThreshold, float errorThreshold2, bool ibd2, bool printCoef, int numThreads, float min_coef, float fudgeFactor){	

	//omp_lock_t lock;//Lock for the coefficient/relationship class file
	//omp_init_lock(&lock);



	float totalGeneticLength = 0.0;//Value to use as whole input genetic length for calculating coefficients.
	for(int chr = 0; chr < Marker::getNumChroms(); chr++)
		totalGeneticLength += Marker::getMarker( Marker::getLastMarkerNum( chr ) )->getMapPos();
	float ibd1Thresholds[9];//Thresholds for class calling
	float ibd2Thresholds[9];
	ibd1Thresholds[0]=0;
	ibd1Thresholds[1]=0.365;//Hardcode for IBD1 only 1st threshold due to not following the formula.
	for(int i = 2; i<=8; i++)
		ibd1Thresholds[i] = 1-(1.0/pow(2.0,((i-1)+1.5)));

	//	{0.365, 1-(1.0/(2.0^(0+1.5))),1-(1.0/(2.0^(1+1.5))),1-(1.0/(2.0^(2+1.5))),1-(1.0/(2.0^(3+1.5))),1-(1.0/(2.0^(4+1.5))),1-(1.0/(2.0^(5+1.5))),1-(1.0/(2.0^(6+1.5))),1-(1.0/(2.0^(7+1.5)))};

	ibd2Thresholds[0]=0.475;
	for(int i = 1; i<=8; i++)
		ibd2Thresholds[i] = (1.0/(pow(2.0,(i+1.5))));

	//Transposes the input matrix and handles the 8->64 bit format conversion..
	interleaveAndConvertData(transposedData, indBlocks, markerWindows, numMarkers, numIndivs, starts, ends,  chromWindowStarts);

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(numThreads); //Forces input specified thread number with default 1. 

	printf("Beginning segment detection with %i thread(s)...",numThreads);
	fflush(stdout);

	/*FileOrGZ<IO_TYPE> pFileArray[numThreads];
	  FileOrGZ<IO_TYPE> classFileArray[numThreads];
	  std::vector<SegmentData> storedSegArray[numThreads];


	  for(int thread = 0; thread < numThreads; thread++){	
	  std::string threadname;
	  threadname = filename +".seg."+std::to_string((1+thread))+extension;
	  bool success = pFileArray[thread].open(threadname.c_str(), "w");
	  if(!success){
	  printf("\nERROR: could not open output VCF file %s!\n", threadname.c_str());
	  perror("open");
	  exit(1);
	  }


	  if(printCoef){
	  std::string classFilename;
	  classFilename = filename +".coef."+std::to_string((1+thread))+extension;
	  success = classFileArray[thread].open(classFilename.c_str(), "w");
	  if(!success){
	  printf("\nERROR: could not open output VCF file %s!\n", classFilename.c_str());
	  perror("open");
	  exit(1);
	  }


	  classFileArray[thread].printf("Individual1\tIndividual2\tKinship_Coefficient\tIBD2_Fraction\tSegment_Count\tDegree\n");
	  }
	  }*/

#pragma omp parallel
	{
		std::string threadname;
		FileOrGZ<IO_TYPE> pFile;
		FileOrGZ<IO_TYPE> classFile;
		threadname = filename +"."+std::to_string((1+omp_get_thread_num()))+".seg"+extension;
		bool success = pFile.open(threadname.c_str(), "w");
		if(!success){
			printf("\nERROR: could not open output VCF file %s!\n", threadname.c_str());
			perror("open");
			exit(1);
		}
		if(printCoef){
			std::string classFilename;
			classFilename = filename +".coef."+std::to_string((1+omp_get_thread_num()))+extension;
			success = classFile.open(classFilename.c_str(), "w");
			if(!success){
				printf("\nERROR: could not open output VCF file %s!\n", classFilename.c_str());
				perror("open");
				exit(1);
			}


			classFile.printf("Individual1\tIndividual2\tKinship_Coefficient\tIBD2_Fraction\tSegment_Count\tDegree\n");
		}

		std::vector<SegmentData> storedSegs;
#pragma omp for schedule(dynamic, 10)
		for(uint64_t indiv1=0; indiv1<numIndivs; indiv1++){
			for(uint64_t indiv2=indiv1+1; indiv2<numIndivs; indiv2++){
				//FileOrGZ<IO_TYPE> pFile = pFileArray[0];
				//FileOrGZ<IO_TYPE> classFile = classFileArray[0];
				//std::vector<SegmentData> storedSegs = storedSegArray[0];
				float totalIBD1Length = 0;
				float totalIBD2Length = 0;

				SegmentData currentIBD1Segment;// = new SegmentData;
				SegmentData currentIBD2Segment;// = new SegmentData;
				currentIBD1Segment.ibdType = 1;


				//IBD2 equivalents of the previously described variables.
				currentIBD2Segment.ibdType = 2;

				HomozygAndMiss *ind1Data = seek(transposedData,indiv1,markerWindows);//Stores the indices of the startpoints for the two individuals in the stored homozygous data and missing data.
				HomozygAndMiss *ind2Data = seek(transposedData,indiv2,markerWindows);


				for(int chr = 0; chr < Marker::getNumChroms(); chr++){
					currentIBD1Segment.chrom = Marker::getChromName(chr);
					currentIBD2Segment.chrom = Marker::getChromName(chr);
					for(uint64_t markerWindow=chromWindowBoundaries[chr][0]; markerWindow<=chromWindowBoundaries[chr][1]; markerWindow++){
						//Booleans for the status of the IBD conclusions for the current window comparison.
						bool isIBD1;
						bool isIBD2;

						int windowErrorCount;
						int windowErrorCount2;


						compareWindowsIBD1(ind1Data->hom1, ind1Data->hom2, ind2Data->hom1, ind2Data->hom2, windowErrorCount, isIBD1);

						//Don't bother checking IBD2 if there are too many errors in IBD1.
						if(ibd2 && isIBD1){
							compareWindowsIBD2(ind1Data->hom1, ind1Data->hom2, ind2Data->hom1, ind2Data->hom2, ind1Data->miss, ind2Data->miss, windowErrorCount2, isIBD2);
						}
						else{
							isIBD2=false;
							windowErrorCount2=3;
						}

						if(ibd2){//The check for if IBD2 is enabled
							if(isIBD2){
								int currentEndPos = ends[markerWindow];
								if(windowErrorCount2>0){
									if(currentIBD2Segment.hasStart() && currentIBD2Segment.errorCheck(windowErrorCount2, currentEndPos, errorThreshold2)){//Update potential endpoint if error threshold passed.
										bool realSeg = endSegment(currentIBD2Segment, min_length2, min_markers2,  storedSegs, totalIBD1Length, totalIBD2Length);
										if(realSeg){
											handleIBD1PostIBD2(currentIBD1Segment,  storedSegs, totalIBD1Length, totalIBD2Length);
										}
										currentIBD2Segment.selfErase();//May be unnecessary relic from earlier draft. Some of these are, some aren't.
									}
									else{//Update endpoints if no errors.
										currentIBD2Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount2, windowErrorCount);
									}

								}
								else{

									if(!currentIBD2Segment.hasStart())//Update startopint if there is no current real segment.
										currentIBD2Segment.updateSegmentStartpoints(starts[markerWindow]);
									currentIBD2Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount2, windowErrorCount);//Update endpoint always.
								}

							}
							else {//If no IBD, attempt to store segment if valid.
								bool realSeg = endSegment(currentIBD2Segment, min_length2, min_markers2,  storedSegs, totalIBD1Length, totalIBD2Length);
								if(realSeg){

									handleIBD1PostIBD2(currentIBD1Segment,  storedSegs, totalIBD1Length, totalIBD2Length);
								}
								currentIBD2Segment.selfErase();
							}
						}




						if(isIBD1){
							int currentEndPos = ends[markerWindow];
							if(windowErrorCount>0){
								if(currentIBD1Segment.hasStart() && currentIBD1Segment.errorCheck(windowErrorCount, currentEndPos, errorThreshold)){
									endSegment(currentIBD1Segment, min_length, min_markers,  storedSegs, totalIBD1Length, totalIBD2Length);
									currentIBD1Segment.selfErase();

								}
								else{
									currentIBD1Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount, windowErrorCount);
								}

							}
							else{
								if(!currentIBD1Segment.hasStart())
									currentIBD1Segment.updateSegmentStartpoints(starts[markerWindow]);
								currentIBD1Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount, windowErrorCount);
							}
						}
						else {
							endSegment(currentIBD1Segment, min_length, min_markers,  storedSegs, totalIBD1Length, totalIBD2Length);
							currentIBD1Segment.selfErase();
						}

						ind1Data++;
						ind2Data++;

					}
					//Handle the end of chromosomes, as there is no end of IBD or error density increase to trigger these segments to print otherwise.

					if(ibd2){
						//Attempt to print IBD2 at the end of the chromosome if there is a valid segment.
						bool realSeg = endSegment(currentIBD2Segment, min_length2, min_markers2,  storedSegs, totalIBD1Length, totalIBD2Length);
						if(realSeg){
							handleIBD1PostIBD2(currentIBD1Segment,  storedSegs, totalIBD1Length, totalIBD2Length);
						}
					}
					//Attempt to print IBD1 at the end of the chromosome if there is a valid segment.
					endSegment(currentIBD1Segment, min_length, min_markers,  storedSegs, totalIBD1Length, totalIBD2Length);
					currentIBD1Segment.selfErase();
					currentIBD2Segment.selfErase();
				}



				//Perform relationship inference and compare to coefficient threshold.
				float relCoef = (0.25*totalIBD1Length+0.5*totalIBD2Length)/totalGeneticLength + fudgeFactor;//NOTE: Currently uses IBD2-based coefficient threshold whether or not IBD2 is enabled.
				float ibd1TotalMod = (1 - (totalIBD1Length + 4 * fudgeFactor)/totalGeneticLength);

				if(relCoef >= min_coef){
					if(printCoef){
						int classVal;
						if(ibd2){
							for(classVal=0; classVal<8; classVal++){
								if((relCoef)>ibd2Thresholds[classVal]){
									break;
								}
							}
						}
						else{
							for(classVal=0; classVal<8; classVal++){
								if(ibd1TotalMod<ibd1Thresholds[classVal]){
									break;
								}
							}
						}
						if(classVal==8)
							classVal=-1;
						//omp_set_lock(&lock);	
						classFile.printf("%s %s %f %f %i %i\n",PersonLoopData::_allIndivs[indiv1]->getId(),PersonLoopData::_allIndivs[indiv2]->getId(), relCoef, storedSegs.size(), totalIBD2Length/totalGeneticLength, classVal);
					}
					//omp_unset_lock(&lock);
					for (SegmentData seg: storedSegs){
						seg.printSegment(pFile,indiv1,indiv2);
					}
				}
				storedSegs.clear();
				//storedSegs.resize(maxSize);
				//delete currentIBD1Segment;
				//delete currentIBD2Segment;

				//internalIteration(transposedData, numIndivs, indBlocks, numMarkers, markerWindows, starts, ends,  filename, extension, chromWindowBoundaries, chromWindowStarts, min_markers, min_length, min_markers2, min_length2, errorThreshold, errorThreshold2, ibd2, printCoef, numThreads, min_coef, fudgeFactor, indiv1, indiv2,  pFileArray, classFileArray, storedSegArray, ibd1Thresholds, ibd2Thresholds, totalGeneticLength);
			}
		}
	}
	//omp_destroy_lock(&lock);
	printf("done.\n");
}




//Consistent Print Statement for usage.
void printUsageAndExit(){

	printf("REQUIRED ARGUMENTS:\n");
	printf(" EITHER:\n\n");
	printf("   First Three Arguments [bed file] [bim file] [fam file]\n");
	printf("      Specifies the plink format files for the data by specific name.\n");
	printf("      Must be first 3 arguments.\n\n");
	printf(" OR:\n\n");
	printf("   -b [prefix] or -bfile [prefix]\n");
	printf("      Specifies the prefix to be used with prefix.bed, prefix.bim, and prefix.fam for the plink format input.\n");
	printf("      Does not need to be first argument.\n\n\n");
	printf("OPTIONS:\n");
	printf(" -mL or -min_l <value>\n");
	printf("      Specify minimum length for acceptible segments to output.\n");
	printf("      Defaults to 7 centimorgans.\n\n");
	printf(" -er or -errorRate <value>\n");
	printf("      Specify acceptible error rate in a segment before considering it false.\n");
	printf("      Defaults to .004 errors per marker.\n\n");
	printf(" -f or file <filename>\n");
	printf("      Specify output file.\n\n");
	printf("      Defaults to ibis<thread number>.seg and will output a separate output file for each thread.\n");
	printf(" -2 or -ibd2\n");
	printf("      Enable ibd2 analyses.\n\n");
	printf(" -mt <value>\n");
	printf("      Set minimum number of markers required for acceptible segments to output.\n");
	printf("      Defaults to 500 markers\n\n");
	printf(" -threads <value>\n");
	printf("      Set the number of threads available to IBIS for parallel processing.\n");
	printf("      Defaults to 1\n\n");
	printf(" -gzip\n");
	printf("      Have the program output gzipped segment files\n\n");
	printf(" -force\n");
	printf("      Prevent the program from trying to convert the input genetic map into centiMorgans\n");
	printf("      Program otherwise multiplies input map by 100 if the input map's total length is less than 6.0 units of genetic distance\n\n");
	printf(" -c\n");
	printf("      Set a minimum kinship coefficient for IBIS to print, omitting pairs of lower relatedness from the output.");
	printf("      Defaults to 0.\n\n");
	printf("OUTPUT FORMAT:\n\n");
	printf(" -Segment file\n");
	printf("sample1 sample2 chrom phys_start_pos phys_end_pos IBD_type gen_start_pos gen_end_pos seg_length marker_count debug_data debug_data\n\n");
	printf(" -Coef file\n");
	printf("sample1 sample2 kinship_coefficient IBD2_fraction segment_count degree_of_relationship\n\n");
	exit(1);


}

int main(int argc, char **argv) {

	const char* VERSION_NUMBER = "1.12";
	const char* RELEASE_DATE = "December 21, 2019";
	printf("IBIS Segment Caller!  v%s    (Released %s)\n\n", VERSION_NUMBER, RELEASE_DATE);

	uint64_t numIndivs, numMarkers;//counts of input quantities.
	float min_length = 7.0, min_length2 = 2.0;//cM minimum lengths
	float min_coef = 0.0;
	bool printCoef = false;
	bool noPrefixGiven = true;//Check if no input file is given.
	bool ibd2 = false;//Check if IBD2 is requested.
	bool bFileNamesGiven = false;//Toggle for input as one argument or 3.
	bool distForce = false;//If true, stops the program from converting the input to cM.
	bool gzip = false;//If True, gzips the output.
	float errorDensityThreshold = 0.004, errorDensityThreshold2 = 0.008;//Maximum allowed error rates.
	float min_markers = 511, min_markers2 = 191;//Marker minimums for segments.
	const char* chrom = NULL;//Checks for specified chromosome. Optional, but required if the input contains multiple chromosomes.
	int numThreads;//Input threadnumber.
	numThreads=0;
	std::string filename, extension;//For building the output file.
	char bfileNameBed[100],bfileNameBim[100],bfileNameFam[100];//locations for input filenames.
	float fudgeFactor = 0.00138;

	printf("Viewing arguments...\n");
	fflush(stdout);
	if(argc == 1){
		printUsageAndExit();		
	}
	for (int i = 0; i < argc; ++i) {
		std::string arg = argv[i];
		if ((arg == "-h") || (arg == "-help")) {
			printUsageAndExit();
		}
		else if( arg == "-b" || arg == "-bfile"){
			bFileNamesGiven = true;
			strcpy(bfileNameBed, argv[i+1]);
			strcat(bfileNameBed, ".bed");
			strcpy(bfileNameBim, argv[i+1]);
			strcat(bfileNameBim, ".bim");
			strcpy(bfileNameFam, argv[i+1]);
			strcat(bfileNameFam, ".fam");
			printf("%s - Running with input files: %s, %s, %s\n",arg.c_str(), bfileNameBed, bfileNameBim, bfileNameFam);

		}
		else if (arg == "-mL" || arg == "-min_l") {
			min_length = atof(argv[i + 1]);
			//min_length2 = min_length/2.0;
			printf("%s - running with minimum IBD1 length %f and minimum IBD2 length %f\n",arg.c_str(),min_length,min_length2);
		}
		else if (arg == "-mL2" || arg == "-min_l2") {
			min_length2 =  atof(argv[i + 1]);
			printf("%s - running with minimum IBD2 length %f\n",arg.c_str(),min_length2);
		}
		else if( arg =="-er" || arg == "-errorRate"){
			errorDensityThreshold=atof(argv[i+1]);
			//errorDensityThreshold2=errorDensityThreshold*2;
			printf("%s - running with error rate %f, IBD2 error rate %f\n",arg.c_str(), errorDensityThreshold, errorDensityThreshold2);
		}
		else if(arg =="-er2" || arg == "-errorRate2"){
			errorDensityThreshold2=atof(argv[i+1]);
			printf("%s - running with error rate %f\n",arg.c_str(), errorDensityThreshold2);

		}
		else if( arg =="-f" || arg == "-file"){
			std::string fileTemp(argv[i+1]);
			filename = fileTemp;
			//pFile = fopen(filename, "w");
			noPrefixGiven = false;
			printf("%s - setting output file %s<thread number>.seg\n",arg.c_str(), filename.c_str());

		}
		else if(arg=="-mt"){
			min_markers=atoi(argv[i+1])-1;
			//min_markers2=min_markers/2.0;
			printf("%s - running with minimum IBD1 marker threshold of %f and minimum IBD2 marker threshold of %f\n",arg.c_str(), min_markers,min_markers2);
		} else if(arg=="-mt2"){
			min_markers2=atoi(argv[i+1])-1;
			printf("%s - running with minimum IBD2 marker threshold of %f\n",arg.c_str(), min_markers2);
		}

		else if(arg =="-chr"){
			chrom = argv[i+1];
			printf("%s - running on chromosome %s\n",arg.c_str(),chrom);
		}
		else if(arg=="-ibd2" || arg == "-2"){
			ibd2=true;
			printf("%s - running with IBD2 detection enabled\n",arg.c_str());
		}
		else if(arg=="-gzip"){
			gzip = true;
			printf("%s - running with gzip enabled\n",arg.c_str());
		}
		else if(arg=="-force"){
			distForce=true;
			printf("%s - forcing morgan format input and output\n",arg.c_str());
		}
		else if(arg=="-threads"){
			numThreads=atoi(argv[i+1]);
			printf("%s - running with %i threads\n",arg.c_str(), numThreads);
		}
		else if(arg=="-c"){
			min_coef = atof(argv[i+1]);
		}
		else if(arg=="-a"){
			fudgeFactor = atoi(argv[i+1]);
		}
		else if(arg=="-printCoef"){
			printCoef = true;
			printf("%s - printing coefficient file\n", arg.c_str());
		}
		else if(arg.at(0)=='-'){
			printUsageAndExit();
		}


	}
	if(noPrefixGiven){
		std::cout<<"No prefix given. Defaulting filenames to \"ibis.<thread number>.seg\"\"\n";
		std::string fileTemp("ibis");
		filename = fileTemp;
	}
	if(bFileNamesGiven)
	{
		PersonIO<PersonLoopData>::readData(bfileNameBed, bfileNameBim, bfileNameFam, /*onlyChr=*/ chrom, /*startPos=*/ 0, /*endPos=*/ INT_MAX, /*XchrName=*/ "X", /*noFamilyId=*/ 1, /*log=*/ NULL, /*allowEmptyParents=*/ false, /*bulkData=*/ false, /*loopData*/ true);
	}
	else{
		printf("No -b or -bfile - Running with input files: %s, %s, %s\n",argv[1], argv[2], argv[3]);
		PersonIO<PersonLoopData>::readData(argv[1], argv[2], argv[3], /*onlyChr=*/ chrom, /*startPos=*/ 0, /*endPos=*/ INT_MAX, /*XchrName=*/ "X", /*noFamilyId=*/ 1, /*log=*/ NULL, /*allowEmptyParents=*/ false, /*bulkData=*/ false, /*loopData*/ true);
	}
	if(numThreads==0){
		numThreads = 1;
	}



	//PersonLoopData::getBulkContainers(dataPointer, bytesPerMarker);
	numIndivs = PersonLoopData::_allIndivs.length();

	//Name change here for error thresholds.
	float errorThreshold = errorDensityThreshold; 
	float errorThreshold2 = errorDensityThreshold2;
	//uint8_t *data = *dataPointer;



	numMarkers = Marker::getNumMarkers();
	uint64_t indBlocks = (numIndivs + 63) / 64; //blocks are based on the number of individuals
	uint64_t markerWindows = 0;

	uint64_t chromBoundaries[Marker::getNumChroms()][2];//Marker indices in the input for the starts and ends of the given chromosomes.
	uint64_t chromWindowBoundaries[Marker::getNumChroms()][2];//Window indices in the stored data for the starts and ends of the given chromosomes.
	uint64_t chromWindowStarts[Marker::getNumChroms()];//Total number of 64-marker blocks that comprise all chromosomes before the indexed chromosome in the stored data. Used to find positions of bits in the data.

	//uint64_t chromMarkerBlocks[Marker::getNumChroms()];

	for (int chrIndex = 0; chrIndex < Marker::getNumChroms(); chrIndex++) {
		chromBoundaries[chrIndex][0]=Marker::getFirstMarkerNum(chrIndex);
		chromBoundaries[chrIndex][1]=Marker::getLastMarkerNum(chrIndex);
	}


	//Lays out total number of blocks of 64 required to represent all markers in the input.
	chromWindowStarts[0] = 0;
	uint64_t markerWindowCount = (Marker::getNumChromMarkers(0) + 63) / 64;
	markerWindows += markerWindowCount;
	for(int chrIndex = 1; chrIndex < Marker::getNumChroms(); chrIndex++){
		chromWindowStarts[chrIndex] = chromWindowStarts[chrIndex-1] + markerWindowCount;
		markerWindowCount = (Marker::getNumChromMarkers(chrIndex) + 63) / 64;
		markerWindows += markerWindowCount;
	}

	if(!distForce){
		for(int chr = 0; chr< Marker::getNumChroms(); chr++){

			if(6.0 > Marker::getMarker( Marker::getLastMarkerNum( chr ) )->getMapPos()){
				printf("Chromosome map shorter than 6 units of genetic distance.\n Morgan input detected - Converting to centimorgans. (Prevent this by running with -force argument)\n");
				Marker::convertMapTocM();
				break;
			}

		}
	}

	float totalGeneticLength = 0.0;//Value to use as whole input genetic length for calculating coefficients.
	for(int chr = 0; chr < Marker::getNumChroms(); chr++)
		totalGeneticLength += Marker::getMarker( Marker::getLastMarkerNum( chr ) )->getMapPos();



	HomozygAndMiss *transposedData = new HomozygAndMiss[numIndivs*markerWindows];




	std::vector<int> starts,ends;


	printf("Defining Windows... ");
	fflush(stdout);
	for(int chr = 0; chr < Marker::getNumChroms(); chr++){
		chromWindowBoundaries[chr][0]=starts.size();
		for(uint64_t x = chromBoundaries[chr][0]; x <= chromBoundaries[chr][1]; x += 64){
			starts.push_back(x);
			ends.push_back(min(chromBoundaries[chr][1], x + 63));
		}
		chromWindowBoundaries[chr][1]=ends.size()-1;
	}


	printf("done.\n");

	if(numThreads>1){
		if(gzip){
			extension = ".gz";
			segmentAnalysis<gzFile>(transposedData, numIndivs, indBlocks, numMarkers, markerWindows, starts, ends, filename, extension, chromWindowBoundaries, chromWindowStarts, min_markers, min_length, min_markers2, min_length2, errorThreshold, errorThreshold2, ibd2, printCoef, numThreads, min_coef, fudgeFactor);	
		}
		else{
			extension = "";
			segmentAnalysis<FILE *>(transposedData, numIndivs, indBlocks, numMarkers, markerWindows, starts, ends, filename, extension, chromWindowBoundaries, chromWindowStarts, min_markers, min_length, min_markers2, min_length2, errorThreshold, errorThreshold2, ibd2, printCoef, numThreads, min_coef, fudgeFactor);	
		}
	}
	else{

		if(gzip){
			extension = ".gz";
			segmentAnalysisMonoThread<gzFile>(transposedData, numIndivs, indBlocks, numMarkers, markerWindows, starts, ends, filename, extension, chromWindowBoundaries, chromWindowStarts, min_markers, min_length, min_markers2, min_length2, errorThreshold, errorThreshold2, ibd2, printCoef, numThreads, min_coef, fudgeFactor);
		}
		else{
			extension = "";
			segmentAnalysisMonoThread<FILE *>(transposedData, numIndivs, indBlocks, numMarkers, markerWindows, starts, ends, filename, extension, chromWindowBoundaries, chromWindowStarts, min_markers, min_length, min_markers2, min_length2, errorThreshold, errorThreshold2, ibd2, printCoef, numThreads, min_coef, fudgeFactor);
		}

	}
}
