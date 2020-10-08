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
#include <assert.h>


//Thresholds for degree calling
//When IBD1 and IBD2 segments called
//These are as in the KING paper: 2^-(degree + 1.5) [except degree 0]
float ibd2Thresholds[9] = { 0.475, 0.176776695, 0.088388348, 0.044194174, 0.022097087, 0.011048543, 0.005524272, 0.002762136 };
//When only IBD status (either IBD1 or 2) is known, how much IBD0? (fractional)
//These are as in the KING paper: degree 1 fixed; others are
// 1 - 2^-(degree + 0.5)
float ibd0Thresholds[9] = { 0, 0.365, 0.6464466, 0.8232233, 0.9116117, 0.9558058, 0.9779029, 0.9889515, 0.9944757 };

template<typename IO_TYPE>
class FileOrGZ {
	public:
		bool open(const char *filename, const char *mode);
		int getline();
		int printf(const char *format, ...);
		size_t pfwrite(const void *ptr, size_t size, size_t n);
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
size_t FileOrGZ<FILE *>::pfwrite(const void *ptr, size_t size, size_t n) {
	return fwrite(ptr, size, n, fp);

}
template<>
size_t FileOrGZ<gzFile>::pfwrite(const void *ptr, size_t size, size_t n) {                                                                                                                                                                      printf("GZip and Binary not currently supported together");
	exit(1);
	return 0;
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
		void addBits(uint64_t h1, uint64_t h2, uint64_t m);//setter function.
		int getMissCount();//getter function
};


void HomozygAndMiss::addBits(uint64_t h1, uint64_t h2, uint64_t m){//Places single-bit values in the input into the bit sets.
	hom1+=h1;
	hom2+=h2;
	miss+=m;
}

int HomozygAndMiss::getMissCount(){//Return count of missing sites for the current sample's window.
	return (int)__builtin_popcountl(miss);
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
		int missingness;//Total count of missing sites in a segment;
		int trackedMissingness;//Missingness beyond the current end of the vlaid segment, extending into the ongoing consecutive error area.
		int ibdType; // Either 1 or 2 depending on nature of IBD for the segment.
		bool realIBD2; // Marks IBD1 segments that were contiguous with a stored real IBD2 segment.
		int chrom; //index of chrom name in the input .bim file.
		static std::vector<float> * altMap;// based on maxDif/maxDist parameter
		SegmentData(){
			startPos = -1;
			startFloat = -1;//genetic start position of ongoing segment.
			errorCount = 0;//error count of ongoing segment.
			trackedError = 0;
			trackedErrorIBD1 = 0;
			missingness = 0;
			trackedMissingness = 0;
			endFloat = -1;
			endPos = -100;//position of last error free window. Used to backtrack segment endings to the last error free window.

			realIBD2 = false;
		}

		template<class IO_TYPE>
			void printSegment(FileOrGZ<IO_TYPE> &pFile, int ind1, int ind2, bool binary);

		template<class IO_TYPE>
			void printHBDSegment(FileOrGZ<IO_TYPE> &pFile, int ind);

		float cMLength();
		int markerLength();
		bool errorCheck(int newError, int newMissingness, int newEndPos, float errorRate);//Confirms segments error
		void updateSegmentEndpoints(int newEndPos, int newError, int newIBD1Error, int newMissingness);
		void updateSegmentStartpoints(int newStartPos);
		bool checkSegment(float min_length, int marker_length);
		bool hasStart();
		void selfErase();
		void trackIBD1Errors();


};

//Declare static member
std::vector<float> * SegmentData::altMap;

inline float SegmentData::cMLength(){return endFloat-startFloat;}
inline int SegmentData::markerLength(){return endPos-startPos+1-missingness;}


template<class IO_TYPE>
void SegmentData::printSegment(FileOrGZ<IO_TYPE> &pFile, int ind1, int ind2, bool binary){
	if(binary){
		int length = markerLength();
		int erC = int(errorCount);
		int startPhys = Marker::getMarker(startPos)->getPhysPos();
		int endPhys = Marker::getMarker(endPos)->getPhysPos();

		//int buffer1[] = {ind1, ind2, chrom, startPhys, endPhys, ibdType, length, erC};

		//PersonLoopData::_allIndivs[ind1]->getId(),PersonLoopData::_allIndivs[ind2]->getId(),Marker::getChromName(chrom), Marker::getMarker(startPos)->getPhysPos(), Marker::getMarker(endPos)->getPhysPos(),ibdType, startFloat, endFloat, cMLength(), markerLength(), int(errorCount), float(errorCount) / float(endPos - startPos)

		//float buffer2[] = {startFloat, endFloat};
		//pFile.pfwrite(buffer1,sizeof(int),8);
		//pFile.pfwrite(buffer2,sizeof(float), 2);
		pFile.pfwrite(&ind1,sizeof(int),1);
		pFile.pfwrite(&ind2,sizeof(int),1);
		pFile.pfwrite(&chrom,sizeof(uint16_t),1);
		pFile.pfwrite(&startPhys,sizeof(int),1);
		pFile.pfwrite(&endPhys,sizeof(int),1);
		pFile.pfwrite(&ibdType, sizeof(uint8_t),1);
		pFile.pfwrite(&startFloat, sizeof(float),1);
		pFile.pfwrite(&endFloat, sizeof(float),1);
		pFile.pfwrite(&length, sizeof(int),1);
		pFile.pfwrite(&erC, sizeof(int),1);
	}
	else
		pFile.printf("%s\t%s\t%s\t%i\t%i\tIBD%i\t%f\t%f\t%f\t%i\t%i\t%f\n", PersonLoopData::_allIndivs[ind1]->getId(),PersonLoopData::_allIndivs[ind2]->getId(),Marker::getChromName(chrom), Marker::getMarker(startPos)->getPhysPos(), Marker::getMarker(endPos)->getPhysPos(),ibdType, startFloat, endFloat, cMLength(), markerLength(), int(errorCount), float(errorCount) / float(markerLength()));

}

template<class IO_TYPE>
void SegmentData::printHBDSegment(FileOrGZ<IO_TYPE> &pFile, int ind) {
	pFile.printf("%s\t%s\t%i\t%i\tHBD\t%f\t%f\t%f\t%i\t%i\t%f\n", PersonLoopData::_allIndivs[ind]->getId(),Marker::getChromName(chrom), Marker::getMarker(startPos)->getPhysPos(), Marker::getMarker(endPos)->getPhysPos(), startFloat, endFloat, cMLength(), markerLength(), int(errorCount), float(errorCount) / float(markerLength()));
}

//Return true if the segment with the new modifications fails the error threshold test: true if fail, false if pass.
inline bool SegmentData::errorCheck(int newError, int newMissingness, int newEndPos, float errorRate){
	return (float(errorCount + trackedError + newError) / float(newEndPos - startPos - missingness - newMissingness)) > errorRate;
}

void SegmentData::updateSegmentEndpoints(int newEndPos, int newError, int newIBD1Error, int newMissingness){
	if(newError<=0){
		endPos = newEndPos;
		endFloat = Marker::getMarker(endPos)->getMapPos();
		errorCount+=trackedError;
		trackedError = 0;
		trackedErrorIBD1 = 0;
		missingness += (trackedMissingness + newMissingness);
		trackedMissingness = 0;
	}
	else{//Won't update real endpoint if there was an error in the window, but tracks the information anyway for error analyses.
		trackedError += newError;
		trackedErrorIBD1 += newIBD1Error;
		trackedMissingness+= newMissingness;
	}
}

void SegmentData::updateSegmentStartpoints(int newStartPos){
	startPos = newStartPos;
	startFloat = Marker::getMarker(startPos)->getMapPos();
	trackedError = 0;
	trackedErrorIBD1 = 0;
	errorCount = 0;
	realIBD2 = false;
	missingness = 0;
	trackedMissingness = 0;
}
//Confirms if segment meets relevant conditions to be a valid segment.
inline bool SegmentData::checkSegment(float min_length, int marker_length){
	return ((endPos > startPos) && (startPos >= 0) && (realIBD2 || (((endPos - startPos - missingness) >= marker_length) && ((altMap->at(endPos) - altMap->at(startPos)) > min_length))));
	//return ((endPos > startPos) && (startPos >= 0) && (realIBD2 || (((endPos - startPos) > marker_length) && ((endFloat - startFloat) > min_length))));
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
	missingness = 0;
	realIBD2 = false;
}


void storeSegment(SegmentData& currentSegment, std::vector<SegmentData> &storedSegments, float &totalIBD1, float &totalIBD2){
	storedSegments.push_back((const SegmentData)(currentSegment));
	totalIBD1+=currentSegment.cMLength() * (currentSegment.ibdType == 1);
	totalIBD2+=currentSegment.cMLength() * (currentSegment.ibdType == 2);
}

void storeHBDSegment(SegmentData& currentSegment, std::vector<SegmentData> &storedSegments, float &totalHBD) {
	storedSegments.push_back((const SegmentData)(currentSegment));
	totalHBD+=currentSegment.cMLength();
}

void storeSegOneBack(SegmentData& currentSegment, std::vector<SegmentData> &storedSegments, float &totalIBD1, float &totalIBD2){
	auto it = storedSegments.end();
	// Insert one behind last segment:
	storedSegments.insert( it - 1, currentSegment );
	totalIBD1+=currentSegment.cMLength() * (currentSegment.ibdType == 1);
	totalIBD2+=currentSegment.cMLength() * (currentSegment.ibdType == 2);
}

//Open the segment and kinship coefficient output files
template<class IO_TYPE>
void openSegCoefOut(FileOrGZ<IO_TYPE> &pFile, FileOrGZ<IO_TYPE> &classFile, FileOrGZ<IO_TYPE> &hbdFile, FileOrGZ<IO_TYPE> &incoefFile, std::string &filename, std::string &extension, bool printCoef, bool binary, bool printFam, bool hbd) {
	std::string threadname;
	std::string binExt;

	bool success;
	if(binary)
		binExt="b";
	else
		binExt="";
	threadname = filename +"."+binExt+"seg"+extension;
	if(binary)
		success = pFile.open(threadname.c_str(), "wb");
	else
		success = pFile.open(threadname.c_str(), "w");

	if(!success){
		printf("\nERROR: could not open output file %s!\n", threadname.c_str());
		perror("open");
		exit(1);
	}
	if(binary){
		uint8_t openVal = 129;

		pFile.pfwrite(&openVal,sizeof(openVal),1);
		int totalChromNameSize = 0;
		for(int chr = 0; chr< Marker::getNumChroms(); chr++){
			totalChromNameSize+=strlen(Marker::getChromName(chr))+1;
		}
		pFile.pfwrite(&totalChromNameSize,sizeof(totalChromNameSize),1);
		for(int chr = 0; chr< Marker::getNumChroms(); chr++){
			pFile.pfwrite(Marker::getChromName(chr),sizeof(char),strlen(Marker::getChromName(chr))+1);
		}

		uint8_t pFam = uint8_t(printFam);
		pFile.pfwrite(&pFam,sizeof(uint8_t),1);
	}

	if(printCoef) {
		std::string classFilename;
		classFilename = filename +".coef"+extension;
		success = classFile.open(classFilename.c_str(), "w");
		if(!success){
			printf("\nERROR: could not open output file %s!\n", classFilename.c_str());
			perror("open");
			exit(1);
		}


		classFile.printf("Individual1\tIndividual2\tKinship_Coefficient\tIBD2_Fraction\tSegment_Count\tDegree\n");
	}

	if (hbd) {
		std::string hbdFilename = filename + ".hbd" + extension;
		success = hbdFile.open(hbdFilename.c_str(), "w");
		if (!success) {
			printf("\nERROR: could not open output file %s!\n", hbdFilename.c_str());
			perror("open");
			exit(1);
		}

		if (printCoef) {
			std::string incoefFilename;
			incoefFilename = filename +".incoef"+extension;
			success = incoefFile.open(incoefFilename.c_str(), "w");
			if(!success) {
				printf("\nERROR: could not open output file %s!\n", incoefFilename.c_str());
				perror("open");
				exit(1);
			}

			incoefFile.printf("Individual\tInbreed_Coefficient\tSegment_Count\n");
		}
	}
}

template<class IO_TYPE>
void closeSegCoefOut(FileOrGZ<IO_TYPE> &pFile, FileOrGZ<IO_TYPE> &classFile, FileOrGZ<IO_TYPE> &hbdFile, FileOrGZ<IO_TYPE> &incoefFile, bool printCoef, bool hbd) {
	pFile.close();
	if (printCoef)
		classFile.close();
	if (hbd) {
		hbdFile.close();
		if (printCoef)
			incoefFile.close();
	}
}

//Segment termination code.
//Tries to merge given segment with the stored "last" segment. If that fails, checks the ongoing segment against all the criteria for segment validity. If it succeeds, it prints the stored segment and stores the ongoing one in the "last" segment slot.
//Returns a boolean describing if the segment in storage is a valid segment according to the thresholds, which, as a side effect, will also always be true if the segment just analyzed was valid.
//template<class IO_TYPE>
bool endSegment(SegmentData& currentSegment, float min_length, int min_markers,  std::vector<SegmentData> &storedSegments, float &totalIBD1, float &totalIBD2){
	bool passedChecks = false;
	if(currentSegment.checkSegment(min_length, min_markers)){
		storeSegment(currentSegment, storedSegments, totalIBD1, totalIBD2);
		passedChecks = true;
	}
	currentSegment.selfErase();
	return passedChecks;


}

//HBD segment termination code.
//Returns a boolean describing if the segment in storage is a valid segment according to the thresholds, which, as a side effect, will also always be true if the segment just analyzed was valid.
bool endHBDSegment(SegmentData& currentSegment, float min_length, int min_markers, std::vector<SegmentData> &storedSegments, float &totalHBD) {
	bool passedChecks = false;
	if(currentSegment.checkSegment(min_length, min_markers)) {
		storeHBDSegment(currentSegment, storedSegments, totalHBD);
		passedChecks = true;
	}
	currentSegment.selfErase();
	return passedChecks;
}

//Breaks and prints underlying IBD1 segment after IBD2 segment confirmed to be real. Also updates IBD1 ongoing segment data.
//template<class IO_TYPE>
//template<class IO_TYPE>
bool handleIBD1PostIBD2(SegmentData &ibd1Segment,  std::vector<SegmentData> &storedSegments, float &totalIBD1, float &totalIBD2, int missingnessCount){
	SegmentData ibd2Segment = storedSegments.back();
	if(ibd2Segment.startPos > ibd1Segment.startPos && ibd1Segment.startPos >= 0){
		ibd1Segment.errorCount = -1;
		int newSegmentTrackedMissingness = ibd1Segment.trackedMissingness;
		int prevSegmentMissingness = (ibd1Segment.missingness + ibd1Segment.trackedMissingness - (ibd2Segment.missingness + ibd2Segment.trackedMissingness - missingnessCount));
		int newSegmentMissingness = ibd1Segment.missingness - prevSegmentMissingness - ibd2Segment.missingness;
		ibd1Segment.updateSegmentEndpoints(ibd2Segment.startPos,0,0,0);
		ibd1Segment.missingness = prevSegmentMissingness;
		int ibd1EndPosHolder = ibd1Segment.endPos;
		storeSegOneBack(ibd1Segment, storedSegments, totalIBD1, totalIBD2);
		ibd1Segment.updateSegmentStartpoints(ibd2Segment.endPos);
		ibd1Segment.updateSegmentEndpoints(ibd1EndPosHolder, 0, 0, 0);
		ibd1Segment.trackedMissingness = newSegmentTrackedMissingness;
		ibd1Segment.errorCount = ibd2Segment.trackedErrorIBD1;//Necessary to put real error count back into ongoing IBD1 segment.
		ibd1Segment.realIBD2 = true;
		ibd1Segment.missingness = newSegmentMissingness;
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

//populate the HBD boolean with true or false depending on if there is valid HBD in the current window and if it had any errors in it.
//Population the errorCount value with the number of heterozygous sites in the window.
void compareWindowHBD(uint64_t hom1, uint64_t hom2, uint64_t missing, int &errorCount, bool &HBD) {
	// want heterozygous sites, which is sites that are _not_ homozygous 1 or 2 or missing
	uint64_t fullCompare = ~(hom1 | hom2 | missing);

	errorCount = __builtin_popcountl(fullCompare);
	HBD = errorCount < 3;
}

uint64_t getMissingness(uint64_t miss1, uint64_t miss2){
	uint64_t fullCompare = (miss1 | miss2);
	return (uint64_t)__builtin_popcountl(fullCompare);
}

//Transpose the input genotypes into a matrix of individuals as rows and blocks of markers as columns to speed up our analyses. Populates transposedData with homozygosity bit sets, and missingData with bit sets of where data is missing.
void interleaveAndConvertData(HomozygAndMiss transposedData[], int64_t indBlocks64, uint64_t markerWindows, uint64_t markers, int64_t individuals, std::vector<int> &starts, std::vector<int> &ends,  uint64_t chromWindowStarts[]) {

	printf("Organizing genotype data for analysis... ");
	fflush(stdout);
	uint64_t dataBlocks = ((individuals+3)/4);//Breakdown of the 8 bit block counts based on the number of individuals in the input. Not the blocks in our own internal data format.

	uint64_t mask1 = 1;

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


				}

			}
		}
	}
	delete[] data;
	PersonIO<PersonLoopData>::closeGeno();
	printf("done.\n");
}

template<class IO_TYPE>
void ibisOn(uint64_t indiv1, std::vector<SegmentData> &storedSegs, HomozygAndMiss transposedData[], int numIndivs, int markerWindows, std::vector<int> &starts, std::vector<int> &ends, uint64_t chromWindowBoundaries[][2], int min_markers, float min_length, int min_markers2, float min_length2, int min_markers_hbd, float min_length_hbd, float errorThreshold, float errorThreshold2, float errorThresholdHBD, bool ibd2, bool hbd, bool printCoef, float min_coef, float fudgeFactor, bool binary, bool printFam, omp_lock_t *lock, FileOrGZ<IO_TYPE> &pFile, FileOrGZ<IO_TYPE> &classFile, FileOrGZ<IO_TYPE> &hbdFile, FileOrGZ<IO_TYPE> &incoefFile, float totalGeneticLength, uint64_t index2End, int &totalSegmentCount) {

	// search for IBD from <indiv1> to all individuals with higher indexes:
	for(uint64_t indiv2 = indiv1+1; indiv2<((uint64_t)numIndivs); indiv2++){
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
			currentIBD1Segment.chrom = chr;
			currentIBD2Segment.chrom = chr;
			for(uint64_t markerWindow=chromWindowBoundaries[chr][0]; markerWindow<=chromWindowBoundaries[chr][1]; markerWindow++){
				//Booleans for the status of the IBD conclusions for the current window comparison.
				bool isIBD1;
				bool isIBD2;

				int windowErrorCount;
				int windowErrorCount2;


				compareWindowsIBD1(ind1Data->hom1, ind1Data->hom2, ind2Data->hom1, ind2Data->hom2, windowErrorCount, isIBD1);
				uint64_t missingnessCount = getMissingness(ind1Data->miss,ind2Data->miss);
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
							if(currentIBD2Segment.hasStart()) {
								if (currentIBD2Segment.errorCheck(windowErrorCount2, missingnessCount, currentEndPos, errorThreshold2)){//Update potential endpoint if error threshold passed.
									bool realSeg = endSegment(currentIBD2Segment, min_length2, min_markers2,  storedSegs, totalIBD1Length, totalIBD2Length);
									if(realSeg){
										handleIBD1PostIBD2(currentIBD1Segment,  storedSegs, totalIBD1Length, totalIBD2Length, missingnessCount);
									}
								}
								else{//Update endpoints if no errors
									currentIBD2Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount2, windowErrorCount, missingnessCount);
								}
							}
						}
						else{

							if(!currentIBD2Segment.hasStart())//Update startopint if there is no current real segment.
								currentIBD2Segment.updateSegmentStartpoints(starts[markerWindow]);
							currentIBD2Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount2, windowErrorCount, missingnessCount);//Update endpoint always.
						}
					}
					else {//If no IBD, attempt to store segment if valid.
						bool realSeg = endSegment(currentIBD2Segment, min_length2, min_markers2,  storedSegs, totalIBD1Length, totalIBD2Length);
						if(realSeg){
							handleIBD1PostIBD2(currentIBD1Segment,  storedSegs, totalIBD1Length, totalIBD2Length, missingnessCount);
						}
					}
				}



				if(isIBD1){
					int currentEndPos = ends[markerWindow];
					if(windowErrorCount>0){
						if(currentIBD1Segment.hasStart()) {
							if (currentIBD1Segment.errorCheck(windowErrorCount, missingnessCount, currentEndPos, errorThreshold)){
								endSegment(currentIBD1Segment, min_length, min_markers,  storedSegs, totalIBD1Length, totalIBD2Length);
							}
							else{
								currentIBD1Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount, windowErrorCount, missingnessCount);
							}
						}
					}
					else{
						if(!currentIBD1Segment.hasStart())
							currentIBD1Segment.updateSegmentStartpoints(starts[markerWindow]);
						currentIBD1Segment.updateSegmentEndpoints(currentEndPos, windowErrorCount, windowErrorCount, missingnessCount);
					}
				}
				else {
					endSegment(currentIBD1Segment, min_length, min_markers,  storedSegs, totalIBD1Length, totalIBD2Length);
				}

				ind1Data++;
				ind2Data++;
			}
			//Handle the end of chromosomes, as there is no end of IBD or error density increase to trigger these segments to print otherwise.

			if(ibd2){
				//Attempt to print IBD2 at the end of the chromosome if there is a valid segment.
				bool realSeg = endSegment(currentIBD2Segment, min_length2, min_markers2,  storedSegs, totalIBD1Length, totalIBD2Length);
				if(realSeg){
					handleIBD1PostIBD2(currentIBD1Segment,  storedSegs, totalIBD1Length, totalIBD2Length,0);
				}
			}
			//Attempt to print IBD1 at the end of the chromosome if there is a valid segment.
			endSegment(currentIBD1Segment, min_length, min_markers,  storedSegs, totalIBD1Length, totalIBD2Length);
		}


		//Perform relationship inference and compare to coefficient threshold.
		float relCoef = (0.25*totalIBD1Length+0.5*totalIBD2Length)/totalGeneticLength + fudgeFactor;//NOTE: Currently uses IBD2-based coefficient threshold whether or not IBD2 is enabled.
		float ibd1TotalMod = (1 - (totalIBD1Length + 4 * fudgeFactor)/totalGeneticLength);

		if(relCoef >= min_coef){//Only prints if current pair is sufficiently related.
			if (lock)
				omp_set_lock(lock);
			if(printCoef){
				int classVal;
				if(ibd2){
					for(classVal=0; classVal<8; classVal++){//find degree of relatedness
						if((relCoef)>ibd2Thresholds[classVal]){
							break;
						}
					}
				}
				else{
					for(classVal=0; classVal<8; classVal++){//find degree of relatedness
						if(ibd1TotalMod<ibd0Thresholds[classVal]){
							break;
						}
					}
				}
				if(classVal==8)
					classVal=-1;//notation for unrelated
				classFile.printf("%s\t%s\t%f\t%f\t%i\t%i\n",PersonLoopData::_allIndivs[indiv1]->getId(),PersonLoopData::_allIndivs[indiv2]->getId(), relCoef, storedSegs.size(), totalIBD2Length/totalGeneticLength, classVal);
			}
			for (SegmentData seg: storedSegs){
				seg.printSegment(pFile,indiv1,indiv2,binary);
				totalSegmentCount+=1;
			}
			if (lock)
				omp_unset_lock(lock);
		}
		storedSegs.clear();
	}

	if (!hbd) {
		// no HBD detection: done
		return;
	}

	// run HBD detection on indiv1
	HomozygAndMiss *ind1Data = seek(transposedData, indiv1, markerWindows);

	float totalHBDLength = 0;
	SegmentData currentHBDSegment;

	for(int chr = 0; chr < Marker::getNumChroms(); chr++) {
		currentHBDSegment.chrom = chr;
		for(uint64_t markerWindow=chromWindowBoundaries[chr][0]; markerWindow<=chromWindowBoundaries[chr][1]; markerWindow++) {

			bool isHBD;
			int windowErrorCount;

			compareWindowHBD(ind1Data->hom1, ind1Data->hom2, ind1Data->miss, windowErrorCount, isHBD);
			uint64_t missingnessCount = (uint64_t)__builtin_popcountl(ind1Data->miss);
			if (isHBD) {
				int currentEndPos = ends[markerWindow];
				if (windowErrorCount > 0) {
					if (currentHBDSegment.hasStart()) {
						if (currentHBDSegment.errorCheck(windowErrorCount, missingnessCount, currentEndPos, errorThresholdHBD))
							endHBDSegment(currentHBDSegment, min_length_hbd, min_markers_hbd, storedSegs, totalHBDLength);
						else
							currentHBDSegment.updateSegmentEndpoints(currentEndPos, windowErrorCount, windowErrorCount, missingnessCount);
					}
				}
				else {
					if(!currentHBDSegment.hasStart())
						currentHBDSegment.updateSegmentStartpoints(starts[markerWindow]);
					currentHBDSegment.updateSegmentEndpoints(currentEndPos, windowErrorCount, windowErrorCount, missingnessCount);
				}
			}
			else {
				endHBDSegment(currentHBDSegment, min_length_hbd, min_markers_hbd, storedSegs, totalHBDLength);
			}

			ind1Data++;
		}
		//Handle the end of chromosomes, as there is no end of IBD or error density increase to trigger these segments to print otherwise.
		endHBDSegment(currentHBDSegment, min_length_hbd, min_markers_hbd, storedSegs, totalHBDLength);
	}

	if (printCoef && (min_coef == 0.0 || totalHBDLength > 0.0)) {
		// Print inbreeding coefficient if either the minimum kinship coefficient is 0 (thus all pairs printed) or if there's some HBD length
		float inbreedCoef = totalHBDLength / totalGeneticLength;
		incoefFile.printf("%s\t%f\t%i\n",PersonLoopData::_allIndivs[indiv1]->getId(), inbreedCoef, storedSegs.size());
	}
	for (SegmentData seg: storedSegs) {
		seg.printHBDSegment(hbdFile,indiv1);
	}
	storedSegs.clear();
}

//Loop over the individuals and windows and perform the segment analysis.
//Used only for monothreaded case to avoid overhead of OMP.
template<class IO_TYPE>
void segmentAnalysisMonoThread(HomozygAndMiss transposedData[], int numIndivs, int indBlocks, int numMarkers, int markerWindows, std::vector<int> &starts, std::vector<int> &ends,  std::string &filename, std::string &extension, uint64_t chromWindowBoundaries[][2], uint64_t chromWindowStarts[], int min_markers, float min_length, int min_markers2, float min_length2, int min_markers_hbd, float min_length_hbd, float errorThreshold, float errorThreshold2, float errorThresholdHBD, bool ibd2, bool hbd, bool printCoef, int numThreads, float min_coef, float fudgeFactor, int index1Start, int index1End, bool binary, float totalGeneticLength, bool printFam){

	//Transposes the input matrix and handles the 8->64 bit format conversion..
	interleaveAndConvertData(transposedData, indBlocks, markerWindows, numMarkers, numIndivs, starts, ends,  chromWindowStarts);
	int totalSegmentCount = 0;

	printf("Beginning segment detection with %i thread(s)...",numThreads);
	fflush(stdout);



	FileOrGZ<IO_TYPE> pFile;
	FileOrGZ<IO_TYPE> classFile;
	FileOrGZ<IO_TYPE> hbdFile;
	FileOrGZ<IO_TYPE> incoefFile;

	openSegCoefOut(pFile, classFile, hbdFile, incoefFile, filename, extension, printCoef, binary, printFam, hbd);

	if(index1End==0)
		index1End=numIndivs-1;
	std::vector<SegmentData> storedSegs;
	for(uint64_t indiv1=index1Start; indiv1<=((uint64_t)index1End); indiv1++){
		ibisOn(indiv1, storedSegs, transposedData, numIndivs, markerWindows, starts, ends, chromWindowBoundaries, min_markers, min_length, min_markers2, min_length2, min_markers_hbd, min_length_hbd, errorThreshold, errorThreshold2, errorThresholdHBD, ibd2, hbd, printCoef, min_coef, fudgeFactor, binary, printFam, /*lock=None=*/ NULL, pFile, classFile, hbdFile, incoefFile, totalGeneticLength, index1End, totalSegmentCount);
	}

	closeSegCoefOut(pFile, classFile, hbdFile, incoefFile, printCoef, hbd);
	printf("Total Segments Found: %i\n",totalSegmentCount);
	printf("done.\n");
}


//Loop over the individuals and windows and perform the segment analysis.
//Includes OMP multithreading
template<class IO_TYPE>
void segmentAnalysis(HomozygAndMiss transposedData[], int numIndivs, int indBlocks, int numMarkers, int markerWindows, std::vector<int> &starts, std::vector<int> &ends,  std::string &filename, std::string &extension, uint64_t chromWindowBoundaries[][2], uint64_t chromWindowStarts[], int min_markers, float min_length, int min_markers2, float min_length2, int min_markers_hbd, float min_length_hbd, float errorThreshold, float errorThreshold2, float errorThresholdHBD, bool ibd2, bool hbd, bool printCoef, int numThreads, float min_coef, float fudgeFactor,int index1Start,int index1End, bool binary, float totalGeneticLength, bool printFam){

	omp_lock_t lock;//Lock for the coefficient/relationship class file
	omp_init_lock(&lock);
	int totalSegmentCount = 0;

	//Transposes the input matrix and handles the 8->64 bit format conversion..
	interleaveAndConvertData(transposedData, indBlocks, markerWindows, numMarkers, numIndivs, starts, ends,  chromWindowStarts);

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(numThreads); //Forces input specified thread number with default 1.

	printf("Beginning segment detection with %i thread(s)...",numThreads);
	fflush(stdout);


	FileOrGZ<IO_TYPE> pFile;
	FileOrGZ<IO_TYPE> classFile;
	FileOrGZ<IO_TYPE> hbdFile;
	FileOrGZ<IO_TYPE> incoefFile;

	openSegCoefOut(pFile, classFile, hbdFile, incoefFile, filename, extension, printCoef, binary, printFam, hbd);

	if(index1End==0)
		index1End = numIndivs;

#pragma omp parallel
	{

		std::vector<SegmentData> storedSegs;
#pragma omp for schedule(dynamic, 60)
		for(uint64_t indiv1 = index1Start; indiv1<((uint64_t)index1End); indiv1++){
			ibisOn(indiv1, storedSegs, transposedData, numIndivs, markerWindows, starts, ends, chromWindowBoundaries, min_markers, min_length, min_markers2, min_length2, min_markers_hbd, min_length_hbd, errorThreshold, errorThreshold2, errorThresholdHBD, ibd2, hbd, printCoef, min_coef, fudgeFactor, binary, printFam, &lock, pFile, classFile, hbdFile, incoefFile, totalGeneticLength, index1End, totalSegmentCount);
		}
	}

	closeSegCoefOut(pFile, classFile, hbdFile, incoefFile, printCoef, hbd);
	printf("Total Segments Found: %i\n",totalSegmentCount);
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
	printf(" Execution options\n");
	printf("  -2 or -ibd2\n");
	printf("      Enable IBD2 analyses.\n");
	printf("  -hbd\n");
	printf("      Enable HBD segment detection.\n");
	printf("  -chr <value>\n");
	printf("      Set specific single chromosome to analyse in an input with multiple chromosomes\n");
	printf("      Defaults to processing all chromosomes in the input\n");
	printf("  -t <value> or -threads <value>\n");
	printf("      Set the number of threads available to IBIS for parallel processing.\n");
	printf("      Defaults to 1\n");
	printf("  -noConvert\n");
	printf("      Prevent IBIS from attempting to convert putative Morgan genetic positions to centiMorgans by multiplying these by 100\n");
	printf("      IBIS makes this conversion if any input chromosome is <= 6 genetic units in length, -noConvert disables\n");
	printf("  -maxDist <value>\n");
	printf("      Set a maximum separation distance between SNPs in the input map.\n");
	printf("      Defaults to being inactive.\n");
	printf("  -setIndexStart <value>\n");
	printf("      Set a start index for the set of samples to be compared against all other samples that apepar later in the input dataset.\n");
	printf("      Must be greater than or equal to 0.\n");
	printf("      Defaults to 0, the index of the first sample.\n");
	printf("  -setIndexEnd <value>\n");
        printf("      Set an end index for the set of samples to be compared against all other samples later than the start index in the input dataset. Includes given index.\n");
	printf("      Must be less than or equal to n-1, where n is the number of samples in the input.\n");
        printf("      Defaults to n-1, the index of the last sample.\n\n");
	printf(" IBD threshold parameters:\n");
	printf("  -er or -errorRate <value>\n");
	printf("      Specify acceptable error rate in a segment before considering it false.\n");
	printf("      Defaults to .004 errors per marker.\n");
	printf("  -mL or -min_l <value>\n");
	printf("      Specify minimum length for acceptable segments to output.\n");
	printf("      Defaults to 7 centimorgans.\n");
	printf("  -mt <value>\n");
	printf("      Set minimum number of markers required for acceptable segments to output.\n");
	printf("      Defaults to 436 markers\n\n");
	printf(" IBD2 threshold parameters: (use with -2 or -ibd2)\n");
	printf("  -er2 or -errorRate2 <value>\n");
	printf("      Specify acceptable error rate in an IBD2 segment before considering it false.\n");
	printf("      Defaults to .008 errors per marker.\n");
	printf("  -mL2 or -min_l2 <value>\n");
	printf("      Specify minimum length for acceptable IBD2 segments to output.\n");
	printf("      Defaults to 2 centimorgans.\n");
	printf("  -mt2 <value>\n");
	printf("      Set minimum number of markers required for acceptable IBD2 segments to output.\n");
	printf("      Defaults to 186 markers\n\n");
	printf(" HBD threshold parameters: (use with -hbd)\n");
	printf("  -erH or -errorRateH <value>\n");
	printf("      Specify acceptable error rate in an HBD segment before considering it false.\n");
	printf("      Defaults to .008 errors per marker.\n");
	printf("  -mLH or -min_lH <value>\n");
	printf("      Specify minimum length for acceptable HBD segments to output.\n");
	printf("      Defaults to 3 centimorgans.\n");
	printf("  -mtH <value>\n");
	printf("      Set minimum number of markers required for acceptable HBD segments to output.\n");
	printf("      Defaults to 186 markers\n\n");
	printf(" Output controls:\n");
	printf("  -f <filename> or -o <filename> or -file <filename>\n");
	printf("      Specify output file prefix.\n");
	printf("      Defaults to ibis\n");
	printf("  -bin or -binary\n");
	printf("      Have the program print the .seg file in binary format. Requires tool to interpret.\n");
	printf("  -gzip\n");
	printf("      Have the program output gzipped segment files\n");
	printf("  -noFamID\n");
	printf("      Have the program omit family IDs from the output, including only individual IDs.\n\n");
	printf(" Kinship and inbreeding coefficient file options:\n");
	printf("  -printCoef\n");
	printf("      Have IBIS print a .coef file and a .incoef (with -hbd) file.\n");
	printf("  -a <value>\n");
	printf("      Set a supplemental factor for IBIS to add to kinship coefficients and use for degree classification.\n");
	printf("      Defaults to 0.00138.\n");
	printf("  -d or -degree <value>\n");
	printf("      Set a minimum degree of relatedness for IBIS to print, omitting pairs of lower relatedness from the output.\n");
	printf("      Defaults to including all degrees.\n");
	printf("      Mutually exclusive with -c\n");
	printf("  -c <value>\n");
	printf("      Set a minimum kinship coefficient for IBIS to print, omitting pairs of lower relatedness from the output.\n");
	printf("      Defaults to 0.\n");
	printf("      Mutually exclusive with -d\n");
	exit(1);


}

int main(int argc, char **argv) {

	const char* VERSION_NUMBER = "1.20.7";
	const char* RELEASE_DATE = "September 17, 2020";
	printf("IBIS Segment Caller!  v%s    (Released %s)\n\n", VERSION_NUMBER, RELEASE_DATE);

	uint64_t numIndivs, numMarkers;//counts of input quantities.
	float min_length = 7.0, min_length2 = 2.0, min_length_hbd = 3.0;//cM minimum lengths
	float min_coef = 0.0;//coefficient threshold for output.
	bool printCoef = false;
	bool noPrefixGiven = true;//Check if no input file is given.
	bool ibd2 = false;//Check if IBD2 is requested.
	bool hbd = false; // detect HBD segments
	bool bFileNamesGiven = false;//Toggle for input as one argument or 3.
	bool distForce = false;//If true, stops the program from converting the input to cM.
	bool gzip = false;//If True, gzips the output.
	float errorThreshold = 0.004, errorThreshold2 = 0.008, errorThresholdHBD = 0.008;//Maximum allowed error rates.
	float min_markers = 435, min_markers2 = 185, min_markers_hbd = 185;//Marker minimums for segments.
	const char* chrom = NULL;//Which chromosome to analyze? If null, analyzes all input chromosomes
	int numThreads;//Input threadnumber.
	numThreads=0;
	bool noFams = 0;
	bool modifiedCoef = false;//Check if -c or -d argument was already used in the input.
	std::string filename, extension;//For building the output file.
	char *bfileNameBed = NULL, *bfileNameBim = NULL, *bfileNameFam = NULL;//locations for input filenames.
	float fudgeFactor = 0.00138;
	int index1Start = 0;
	int index1End = 0;
	printf("Viewing arguments...\n");
	fflush(stdout);
	bool binary = false;
	float maxDif = 100000;
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
			uint32_t filenameLength = strlen(argv[i+1]) + 4 + 1; // +1 for '\0'
			bfileNameBed = new char[filenameLength];
			bfileNameBim = new char[filenameLength];
			bfileNameFam = new char[filenameLength];
			sprintf(bfileNameBed, "%s.bed", argv[i+1]);
			sprintf(bfileNameBim, "%s.bim", argv[i+1]);
			sprintf(bfileNameFam, "%s.fam", argv[i+1]);
			printf("%s - Running with input files: %s, %s, %s\n",arg.c_str(), bfileNameBed, bfileNameBim, bfileNameFam);

		}
		else if (arg == "-mL" || arg == "-min_l") {
			min_length = atof(argv[i + 1]);
			printf("%s - running with minimum IBD1 length %f\n",arg.c_str(),min_length);
		}
		else if (arg == "-mL2" || arg == "-min_l2") {
			min_length2 =  atof(argv[i + 1]);
			printf("%s - running with minimum IBD2 length %f\n",arg.c_str(),min_length2);
		}
		else if (arg == "-mLH" || arg == "-min_lH") {
			min_length_hbd = atof(argv[i + 1]);
			printf("%s - running with minimum HBD length %f\n",arg.c_str(),min_length_hbd);
		}
		else if( arg =="-er" || arg == "-errorRate"){
			errorThreshold=atof(argv[i+1]);
			printf("%s - running with error rate %f\n",arg.c_str(), errorThreshold);
		}
		else if(arg =="-er2" || arg == "-errorRate2"){
			errorThreshold2=atof(argv[i+1]);
			printf("%s - running with IBD2 error rate %f\n",arg.c_str(), errorThreshold2);

		}
		else if (arg == "-erH" || arg == "-errorRateHBD") {
			errorThresholdHBD = atof(argv[i+1]);
			printf("%s - running with HBD error rate %f\n",arg.c_str(), errorThresholdHBD);
		}
		else if( arg =="-f" || arg == "-file" || arg == "-o"){
			std::string fileTemp(argv[i+1]);
			filename = fileTemp;
			noPrefixGiven = false;
			printf("%s - setting output file %s.seg\n",arg.c_str(), filename.c_str());

		}
		else if(arg=="-mt"){
			min_markers=atoi(argv[i+1])-1;
			printf("%s - running with minimum IBD1 marker threshold of %f\n",arg.c_str(), min_markers);
		} else if(arg=="-mt2"){
			min_markers2=atoi(argv[i+1])-1;
			printf("%s - running with minimum IBD2 marker threshold of %f\n",arg.c_str(), min_markers2);
		}
		else if (arg == "-mtH"){
			min_markers_hbd = atoi(argv[i+1])-1;
			printf("%s - running with minimum HBD marker threshold of %f\n",arg.c_str(), min_markers_hbd);
		}

		else if(arg =="-chr"){
			chrom = argv[i+1];
			printf("%s - running on chromosome %s\n",arg.c_str(),chrom);
		}
		else if(arg=="-ibd2" || arg == "-2"){
			ibd2=true;
			printf("%s - running with IBD2 detection enabled\n",arg.c_str());
		}
		else if (arg == "-hbd") {
			hbd = true;
			printf("%s - runnning with HBD detection enabled\n",arg.c_str());
		}
		else if(arg=="-gzip"){
			gzip = true;
			printf("%s - running with gzip enabled\n",arg.c_str());
		}
		else if(arg=="-noConvert"){
			distForce=true;
			printf("%s - forcing morgan format input and output\n",arg.c_str());
		}
		else if(arg=="-threads" || arg== "-t"){
			numThreads=atoi(argv[i+1]);
			printf("%s - running with %i threads\n",arg.c_str(), numThreads);
		}
		else if(arg=="-c"){
			if(modifiedCoef){
				printf("Cannot use both -c and -d, or use either more than once\n");
				exit(1);
			}
			modifiedCoef=true;
			min_coef = atof(argv[i+1]);
		}
		else if(arg=="-a"){
			fudgeFactor = atoi(argv[i+1]);
		}
		else if(arg=="-bin"||arg=="-binary"){
			printf("%s - setting binary output.\n",arg.c_str());
			binary=true;
		}
		else if(arg=="-printCoef"){
			printCoef = true;
			printf("%s - printing coefficient file\n", arg.c_str());
		}
		else if(arg=="-noFamID"){
			noFams=1;
			printf("%s - assuming no Family ID in input\n", arg.c_str());
		}
		else if(arg=="-degree" || arg=="-d"){
			if(modifiedCoef){
                                printf("Cannot use both -c and -d, or use either more than once\n");
                                exit(1);
                        }
                        modifiedCoef=true;
			min_coef = (1.0/(pow(2.0,(atoi(argv[i+1])+1.5))));
			printf("%s - creating min kinship coefficient threshold from degree %s as %f\n", arg.c_str(),argv[i+1],min_coef);
		}
		else if(arg=="-setIndexEnd"){
			index1End = atoi(argv[i+1]);
			printf("%s - attempting to set a sample end index of %i\n",arg.c_str(),index1End);
			//Does not yet know size of sample set. Checks later.
			if(index1End<=0){
				printf("Cannot select starting ending index less than or equal to 0.");
				exit(1);
			}
		}
		else if(arg=="-setIndexStart"){
			index1Start = atoi(argv[i+1]);
			printf("%s - attempting to set a sample start index of %i\n",arg.c_str(),index1Start);
			if(index1Start < 0){
				printf("Cannot select starting index for samples less than 0.");
				exit(1);
			}
		}
		else if(arg=="-maxDist"){
			printf("%s - altering maximum SNP separation value to %s\n", arg.c_str(),argv[i+1]);
			maxDif = atof(argv[i+1]);
		}
		else if(arg.at(0)=='-'){
			printf("Unrecognized Argument: %s\n",arg.c_str());
			printUsageAndExit();
		}


	}
	if(noPrefixGiven){
		std::cout<<"No prefix given. Defaulting filenames to prefix \"ibis\"\n";
		std::string fileTemp("ibis");
		filename = fileTemp;
	}
	if(bFileNamesGiven)
	{
		PersonIO<PersonLoopData>::readData(bfileNameBed, bfileNameBim, bfileNameFam, /*onlyChr=*/ chrom, /*startPos=*/ 0, /*endPos=*/ INT_MAX, /*XchrName=*/ "X", /*noFamilyId=*/ noFams, /*log=*/ NULL, /*allowEmptyParents=*/ false, /*bulkData=*/ false, /*loopData*/ true,  /*useParents=*/ false, /*ignoreAlleles=*/ true);
	}
	else{
		printf("No -b or -bfile - Running with input files: %s, %s, %s\n",argv[1], argv[2], argv[3]);
		PersonIO<PersonLoopData>::readData(argv[1], argv[2], argv[3], /*onlyChr=*/ chrom, /*startPos=*/ 0, /*endPos=*/ INT_MAX, /*XchrName=*/ "X", /*noFamilyId=*/ noFams, /*log=*/ NULL, /*allowEmptyParents=*/ false, /*bulkData=*/ false, /*loopData*/ true, /*useParents=*/ false, /*ignoreAlleles=*/ true);
	}



	if(numThreads==0){
		numThreads = 1;
	}

	if(Marker::getNumChroms() > UINT16_MAX-1){
		printf("IBIS cannot handle more than %d different chromosomes in the input due to formatting issues. Please separate them and run ibis on subsets of these.\n",UINT16_MAX-1);
		exit(1);
	}

	numIndivs = PersonLoopData::_allIndivs.length();
	if((uint64_t)index1End > numIndivs-1){
                printf("Cannot select ending index for samples beyond the final index.\n");
        	exit(1);
        }
	else if ((index1End < index1Start) && index1End!=0){
		printf("Cannot select ending index before starting index.\n");
		exit(1);
	}
	numMarkers = Marker::getNumMarkers();
	uint64_t indBlocks = (numIndivs + 63) / 64; //blocks are based on the number of individuals
	uint64_t markerWindows = 0;

	uint64_t chromWindowBoundaries[Marker::getNumChroms()][2];//Window indices in the stored data for the starts and ends of the given chromosomes.
	uint64_t chromWindowStarts[Marker::getNumChroms()];//Total number of 64-marker blocks that comprise all chromosomes before the indexed chromosome in the stored data. Used to find positions of bits in the data.
	SegmentData::altMap = new std::vector<float>;

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
				printf("Chromosome map shorter than 6 units of genetic distance.\n Morgan input detected - Converting to centimorgans. (Prevent this by running with -noConvert argument)\n");
				Marker::convertMapTocM();
				break;
			}

		}
	}

	SegmentData::altMap->push_back(Marker::getMarker(0)->getMapPos());
	//Set up alternative map for max distance control.
	for(int altmark = 1; altmark < Marker::getNumMarkers(); altmark++){

		if(Marker::getMarker(altmark-1)->getChromIdx()==Marker::getMarker(altmark)->getChromIdx()){
			float diff = Marker::getMarker(altmark)->getMapPos()-Marker::getMarker(altmark-1)->getMapPos();
			float mv1 = SegmentData::altMap->back()+maxDif;
			float mv2 = SegmentData::altMap->back()+diff;
			if(mv1>mv2)
				SegmentData::altMap->push_back(mv2);
			else
				SegmentData::altMap->push_back(mv1);

		}
		else
			SegmentData::altMap->push_back(Marker::getMarker(altmark)->getMapPos());
	}

	HomozygAndMiss *transposedData = new HomozygAndMiss[numIndivs*markerWindows];//Represents the genotypes in the format to be used by IBIS.




	std::vector<int> starts,ends;

	//Delineate window boundaries for bitwise comparisons.
	printf("Defining Windows... ");
	fflush(stdout);
	for(int chr = 0; chr < Marker::getNumChroms(); chr++){
		chromWindowBoundaries[chr][0]=starts.size();
		for(uint64_t x = Marker::getFirstMarkerNum(chr); x <= (uint64_t) Marker::getLastMarkerNum(chr); x += 64){
			starts.push_back(x);
			ends.push_back(min(Marker::getLastMarkerNum(chr), x + 63));
		}
		chromWindowBoundaries[chr][1]=ends.size()-1;
	}


	printf("done.\n");

	float totalGeneticLength = 0.0;//Value to use as whole input genetic length for calculating coefficients.
	for(int chr = 0; chr < Marker::getNumChroms(); chr++)
		totalGeneticLength += (Marker::getMarker( Marker::getLastMarkerNum( chr ) )->getMapPos() - Marker::getMarker( Marker::getFirstMarkerNum( chr ) )->getMapPos());
	printf("Total Genetic Length in use:%f\n", totalGeneticLength);
	if(numThreads>1){
		if(gzip){
			extension = ".gz";
			segmentAnalysis<gzFile>(transposedData, numIndivs, indBlocks, numMarkers, markerWindows, starts, ends, filename, extension, chromWindowBoundaries, chromWindowStarts, min_markers, min_length, min_markers2, min_length2, min_markers_hbd, min_length_hbd, errorThreshold, errorThreshold2, errorThresholdHBD, ibd2, hbd, printCoef, numThreads, min_coef, fudgeFactor, index1Start,index1End, binary, totalGeneticLength, noFams);
		}
		else{
			extension = "";
			segmentAnalysis<FILE *>(transposedData, numIndivs, indBlocks, numMarkers, markerWindows, starts, ends, filename, extension, chromWindowBoundaries, chromWindowStarts, min_markers, min_length, min_markers2, min_length2, min_markers_hbd, min_length_hbd, errorThreshold, errorThreshold2, errorThresholdHBD, ibd2, hbd, printCoef, numThreads, min_coef, fudgeFactor, index1Start,index1End, binary, totalGeneticLength, noFams);
		}
	}
	else{
		//Monothreaded to avoid potential overhead. Bad style but a runtime gain.
		if(gzip){
			extension = ".gz";
			segmentAnalysisMonoThread<gzFile>(transposedData, numIndivs, indBlocks, numMarkers, markerWindows, starts, ends, filename, extension, chromWindowBoundaries, chromWindowStarts, min_markers, min_length, min_markers2, min_length2, min_markers_hbd, min_length_hbd, errorThreshold, errorThreshold2, errorThresholdHBD, ibd2, hbd, printCoef, numThreads, min_coef, fudgeFactor,index1Start,index1End, binary, totalGeneticLength, noFams);
		}
		else{
			extension = "";
			segmentAnalysisMonoThread<FILE *>(transposedData, numIndivs, indBlocks, numMarkers, markerWindows, starts, ends, filename, extension, chromWindowBoundaries, chromWindowStarts, min_markers, min_length, min_markers2, min_length2, min_markers_hbd, min_length_hbd, errorThreshold, errorThreshold2, errorThresholdHBD, ibd2, hbd, printCoef, numThreads, min_coef, fudgeFactor,index1Start,index1End, binary, totalGeneticLength, noFams);
		}

	}
}
