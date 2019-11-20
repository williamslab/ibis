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



//Find the lowest order non-zero bit in the uint64_t, and return it in a uint64_t with only that bit set to 1.
inline uint64_t lowestSetBit(uint64_t bitSet)
{
	return (-bitSet) & bitSet;
};

//Determine if the two given segments overlap
inline bool segmentOverlap(int start1, int start2, int end1, int end2){
	return (end2 >= start1) && (end1 >= start2);//(start1 <= end1 && start1 >= start2) || (end1<=end2 && end1>=start2);

}


//Future cosmetic project for print statements. May want later.
/*void sendToStream(ogztream *outStream, char *ind1, char* ind2, int ibdType, const char* chrom, int physStart,int physEnd, char* snpStart, char *snpEnd, float genStart, float genEnd, float length, int markerLength, int lastErrorCount, float errorDensity){
  outStream << ind1<< "\t"<< ind2<< "\t"<< ibdType<< "\t"<< chrom<< "\t"<< physStart<< "\t"<< physEnd<< "\t"<< snpStart<< "\t"<< snpEnd<< "\t"<< genStart<< "\t"<< genEnd< "\t"<< length<< "\t"<< markerLength<< "\t"<< lastErrorCount<< "\t"<< errorDensity<<"\n";
  }*/


//Attempts to merge the segment in storage with the ongoing segment, returns true if able. Returns false if either unable or if merge is not enabled.
//The rest of the code is designed to not print a segment if a valid merge occurs, as future merges may be necessary, and a lot of processing is only triggered by failure to merge.
bool mergeSegments(int &startPos, int &endPos, int &lastStartPos, int &lastEndPos, float &startFloat, float &endFloat, float &lastStartFloat, float &lastEndFloat, float min_length, float min_markers, float &errorCount, float &lastErrorCount, int &lastWindowErrorCount, int ibdType, bool &realIBD2, bool &realIBD2New, bool noMerge){

	if(noMerge)
		return false;//Always the result at the moment.
	float min_l = 2.0;//TODO Make this a function of the input?
	if ((startFloat!=-1) && (startPos-lastEndPos < 66) && (lastWindowErrorCount < 3) && ((ibdType==1 && realIBD2) || ((lastEndFloat-lastStartFloat)>(min_l))) && ((ibdType==1 && realIBD2New) || (endFloat-startFloat >(min_l)))){
		lastEndFloat=endFloat;
		lastEndPos=endPos;
		lastErrorCount+=errorCount+lastWindowErrorCount;
		realIBD2 = realIBD2 || realIBD2New;
		realIBD2New = false;
		return true;
	}
	return false;
}



//Segment termination code.
//Tries to merge given segment with the stored "last" segment. If that fails, checks the ongoing segment against all the criteria for segment validity. If it succeeds, it prints the stored segment and stores the ongoing one in the "last" segment slot.
//Returns a boolean describing if the segment in storage is a valid segment according to the thresholds, which, as a side effect, will also always be true if the segment just analyzed was valid.
template<class IO_TYPE>
bool endSegment(int &startPos, int &endPos, int &lastStartPos, int &lastEndPos, float &startFloat, float &endFloat, float &lastStartFloat, float &lastEndFloat, float min_length, float min_markers, float &errorCount, float &lastErrorCount, int &lastWindowErrorCount, int ibdType, const char* chrom, uint64_t indiv1, uint64_t indiv2, std::vector<std::string> &names, std::vector<std::string> &SNPID, std::vector<float> &geneMap, std::vector<int> &physMap, FileOrGZ<IO_TYPE> &pFile, bool &realIBD2, bool &realIBD2New, bool noMerge){
	bool validMerge = false;
	if(!mergeSegments(startPos, endPos, lastStartPos, lastEndPos, startFloat, endFloat, lastStartFloat, lastEndFloat, min_length, min_markers, errorCount, lastErrorCount, lastWindowErrorCount, ibdType, realIBD2, realIBD2New, noMerge)){
		int segMarkerLength = lastEndPos-lastStartPos;
		float segCMLength = lastEndFloat-lastStartFloat;


		if(lastStartPos !=-100 && segMarkerLength>0 && ((realIBD2 && ibdType==1) || ((segCMLength > min_length) && (segMarkerLength > min_markers)))){
			pFile.printf("%s\t%s\t%s\t%i\t%i\tIBD%i\t%f\t%f\t%f\t%i\t%i\t%f\n", PersonLoopData::_allIndivs[indiv1]->getId(),PersonLoopData::_allIndivs[indiv2]->getId(),chrom, physMap[lastStartPos], physMap[lastEndPos],ibdType, lastStartFloat, lastEndFloat, segCMLength, segMarkerLength, int(lastErrorCount), lastErrorCount / float(segMarkerLength));

		}	

		lastStartPos=startPos;
		lastEndPos=endPos;
		lastStartFloat=startFloat;
		lastEndFloat=endFloat;
		lastErrorCount=errorCount;
		if(ibdType==1){
			realIBD2=realIBD2New;
			realIBD2New=false;
		}


	}
	int segMarkerLength = lastEndPos-lastStartPos;
	float segCMLength = lastEndFloat-lastStartFloat;


	if(lastStartPos!=-100 && (segCMLength > min_length) && (segMarkerLength > min_markers))
		validMerge=true;
	startPos=-100;
	startFloat=-1;
	if(ibdType==1)//If an IBD1 segment was printed, the realIBD2New quantity has been resolved and there is no linked real IBD2 with the ongoing segment from here on out.
		realIBD2New=false;
	errorCount=0;
	return validMerge;

}

//Breaks and prints underlying IBD1 segment after IBD2 segment confirmed to be real. Also updates IBD1 ongoing segment data.
template<class IO_TYPE>
void handleIBD1PostIBD2(int &startPos, int &endPos, int &lastStartPos, int &lastEndPos, float &startFloat, float &endFloat, float &lastStartFloat, float &lastEndFloat, float min_length, float min_markers, float &errorCount, float &lastErrorCount, int &lastWindowErrorCount, int &startPos2, int &endPos2, int &lastStartPos2, int &lastEndPos2, float &startFloat2, float &endFloat2, float &lastStartFloat2, float &lastEndFloat2, float min_length2, float min_markers2, float &errorCount2, float &lastErrorCount2, int &lastWindowErrorCount2, const char* chrom, uint64_t indiv1, uint64_t indiv2, std::vector<std::string> &names, std::vector<std::string> &SNPID, std::vector<float> &geneMap, std::vector<int> &physMap, FileOrGZ<IO_TYPE> &pFile, bool &realIBD2, bool &realIBD2New, bool noMerge){

	if(segmentOverlap(lastStartPos,lastStartPos2,lastEndPos,lastEndPos2)){

		//Checking if a true IBD1 segment needs to be printed before the current IBD2 segment in the stored segment slot is printed.
		if(lastStartPos2>lastStartPos && lastStartPos!=-100){
			int tempStartPos = lastStartPos;
			int tempEndPos = lastStartPos2;
			int segMarkerLength = tempEndPos-tempStartPos;
			float tempStartFloat = geneMap[tempStartPos];
			float tempEndFloat = geneMap[tempEndPos];
			float segCMLength = tempEndFloat-tempStartFloat;
			pFile.printf("%s\t%s\t%s\t%i\t%i\tIBD%i\t%f\t%f\t%f\t%i\t%i\t%f\n", PersonLoopData::_allIndivs[indiv1]->getId(), PersonLoopData::_allIndivs[indiv2]->getId(), chrom, physMap[tempStartPos], physMap[tempEndPos], 1, tempStartFloat, tempEndFloat, segCMLength, segMarkerLength, int(lastErrorCount), lastErrorCount / float(segMarkerLength));


		}
		if(lastEndPos2<lastEndPos && lastEndPos!=-100){
			lastStartPos=lastEndPos2;
			lastStartFloat = lastEndFloat2;
			realIBD2 = true;
		}
		else{
			lastStartPos=-100;
			lastEndPos = -100;
			realIBD2 = false;
		}

	}
	if(segmentOverlap(startPos,lastStartPos2, endPos, lastEndPos2)){
		//checking if a true IBD1 segment needs to be printed before the current IBD2 segment is placed in the stored segment.
		if(startPos!=-100 && lastStartPos2>=startPos){
			int tempStartPos = startPos;
			int tempEndPos = lastStartPos2;
			int segMarkerLength = tempEndPos-tempStartPos;
			float tempStartFloat = geneMap[tempStartPos];
			float tempEndFloat = geneMap[tempEndPos];
			float segCMLength = tempEndFloat-tempStartFloat;
			float segCMLengthOld = lastEndFloat-lastStartFloat;
			int segMarkerLengthOld = lastEndPos-lastStartPos;
			if(!mergeSegments(tempStartPos, tempEndPos, lastStartPos, lastEndPos, tempStartFloat, tempEndFloat, lastStartFloat, lastEndFloat, min_length, min_markers, errorCount, lastErrorCount, lastWindowErrorCount, 1, realIBD2, realIBD2New, noMerge)){
				//Handle possible existing IBD unaffiliated with the IBD2 segment
				if(lastStartPos !=-100 && ((segCMLengthOld > min_length) && (segMarkerLengthOld > min_markers))){
					pFile.printf("%s\t%s\t%s\t%i\t%i\tIBD%i\t%f\t%f\t%f\t%i\t%i\t%f\n", PersonLoopData::_allIndivs[indiv1]->getId(), PersonLoopData::_allIndivs[indiv2]->getId(), chrom, physMap[lastStartPos], physMap[lastEndPos], 1, lastStartFloat, lastEndFloat, segCMLengthOld, segMarkerLengthOld, int(lastErrorCount), lastErrorCount / float(segMarkerLengthOld));
				}
				if(segCMLength>0){//TODO IBD1 error not tracked properly here!
					pFile.printf("%s\t%s\t%s\t%i\t%i\tIBD%i\t%f\t%f\t%f\t%i\t%i\t%f\n", PersonLoopData::_allIndivs[indiv1]->getId(), PersonLoopData::_allIndivs[indiv2]->getId(), chrom, physMap[tempStartPos], physMap[tempEndPos], 1, tempStartFloat, tempEndFloat, segCMLength, segMarkerLength, int(lastErrorCount), lastErrorCount / float(segMarkerLength));
				}
			}
			else{
				segCMLength = lastEndFloat-lastStartFloat;
				segMarkerLength = lastEndPos-lastStartPos;
				pFile.printf("%s\t%s\t%s\t%i\t%i\tIBD%i\t%f\t%f\t%f\t%i\t%i\t%f\n", PersonLoopData::_allIndivs[indiv1]->getId(), PersonLoopData::_allIndivs[indiv2]->getId(), chrom, physMap[lastStartPos], physMap[lastEndPos], 1, lastStartFloat, lastEndFloat, segCMLength, segMarkerLength, int(lastErrorCount), lastErrorCount / float(segMarkerLength));				

			}
			lastStartPos=-100;
			lastEndPos=-100;
			lastStartFloat=-1;
			lastEndFloat=-1;
			realIBD2=false;

		}
		if(startPos!=-100 && lastEndPos2<=endPos){
			startPos=lastEndPos2;
			startFloat = lastEndFloat2;
			realIBD2New=true;//Track the fact that the ongoing IBD1 segment is now linked to true IBD2 and will need to be printed regardless of length later on.
		}
		else{
			startPos=-100;
			startFloat=-1;
			endPos = -100;
			endFloat=-1;
			realIBD2New = false;
		}
	}


}


//CHeck for a valid IBD1 segment in either the stored segment or the ongoing one and print them if valid.
template<class IO_TYPE>
bool forceEndSegment(int &startPos, int &endPos, int &lastStartPos, int &lastEndPos, float &startFloat, float &endFloat, float &lastStartFloat, float &lastEndFloat, float min_length, float min_markers, float &errorCount, float &lastErrorCount, int &lastWindowErrorCount, int ibdType, const char* chrom, uint64_t indiv1, uint64_t indiv2, std::vector<std::string> &names, std::vector<std::string> &SNPID, std::vector<float> &geneMap, std::vector<int> &physMap, FileOrGZ<IO_TYPE> &pFile, bool &realIBD2, bool &realIBD2New, bool noMerge){




	bool printed = false;
	//make one last attempt at merging the ongoing segment with the last one.
	bool merged = mergeSegments(startPos, endPos, lastStartPos, lastEndPos, startFloat, endFloat, lastStartFloat, lastEndFloat, min_length, min_markers, errorCount, lastErrorCount, lastWindowErrorCount, ibdType, realIBD2, realIBD2New, noMerge);
	int segMarkerLength = lastEndPos-lastStartPos;
	float segCMLength = lastEndFloat-lastStartFloat;
	//check the stored segment, print if valid.
	if(lastStartPos!=-100 && segMarkerLength>0 && (realIBD2 || ((segCMLength > min_length) && (segMarkerLength > min_markers)))){
		pFile.printf("%s\t%s\t%s\t%i\t%i\tIBD%i\t%f\t%f\t%f\t%i\t%i\t%f\n", PersonLoopData::_allIndivs[indiv1]->getId(), PersonLoopData::_allIndivs[indiv2]->getId(), chrom, physMap[lastStartPos], physMap[lastEndPos], ibdType, lastStartFloat, lastEndFloat, segCMLength, segMarkerLength, int(lastErrorCount), lastErrorCount / float(segMarkerLength));
		printed=true;
	}
	if((!merged) && (startFloat!=-1)){//If there still is an ongoing segment after merging, check nad print it.
		segMarkerLength = endPos-startPos;
		segCMLength = endFloat-startFloat;
		if(segMarkerLength>0 && ((realIBD2New) || ((segCMLength > min_length) && (segMarkerLength > min_markers)))){
			pFile.printf("%s\t%s\t%s\t%i\t%i\tIBD%i\t%f\t%f\t%f\t%i\t%i\t%f\n", PersonLoopData::_allIndivs[indiv1]->getId(), PersonLoopData::_allIndivs[indiv2]->getId(), chrom, physMap[startPos], physMap[endPos], ibdType, startFloat, endFloat, segCMLength, segMarkerLength, int(errorCount), errorCount / float(segMarkerLength));
			printed=true;
		}
	}
	return printed;
}

//CHeck for a valid IBD2 segment in either the stored segment or the ongoing one and print them if valid. Has additional processing required between printing the IBD2 segment and the associated IBD1 segments.
template<class IO_TYPE>
bool forceEndSegment2(int &startPos2, int &endPos2, int &lastStartPos2, int &lastEndPos2, float &startFloat2, float &endFloat2, float &lastStartFloat2, float &lastEndFloat2, float min_length2, int min_markers2, float &errorCount2, float &lastErrorCount2, int &lastWindowErrorCount2, const char* chrom, uint64_t indiv1, uint64_t indiv2, std::vector<std::string> &names, std::vector<std::string> &SNPID, std::vector<float> &geneMap, std::vector<int> &physMap, FileOrGZ<IO_TYPE> &pFile, bool &realIBD2, bool &realIBD2New, int &startPos, int &endPos, int &lastStartPos, int &lastEndPos, float &startFloat, float &endFloat, float &lastStartFloat, float &lastEndFloat, float min_length, int min_markers, float &errorCount, float &lastErrorCount, int &lastWindowErrorCount,bool noMerge){

	bool printed = false;
	 //make one last attempt at merging the ongoing segment with the last one.
	bool merged = mergeSegments(startPos2, endPos2, lastStartPos2, lastEndPos2, startFloat2, endFloat2, lastStartFloat2, lastEndFloat2, min_length2, min_markers2, errorCount2, lastErrorCount2, lastWindowErrorCount2, 2, realIBD2, realIBD2New, noMerge);
	int segMarkerLength = lastEndPos2-lastStartPos2;
	float segCMLength = lastEndFloat2-lastStartFloat2;
	if(lastStartPos2!=-100 && ((segCMLength > min_length2) && (segMarkerLength > min_markers2))){
		pFile.printf("%s\t%s\t%s\t%i\t%i\tIBD%i\t%f\t%f\t%f\t%i\t%i\t%f\n", PersonLoopData::_allIndivs[indiv1]->getId(), PersonLoopData::_allIndivs[indiv2]->getId(), chrom, physMap[lastStartPos2], physMap[lastEndPos2], 2, lastStartFloat2, lastEndFloat2, segCMLength, segMarkerLength, int(lastErrorCount2), lastErrorCount2 / float(segMarkerLength));
		printed=true;
	}

	if(printed){
		//Deal with the IBD1 contiguous with the processed IBD2 segment.
		handleIBD1PostIBD2<IO_TYPE>(startPos, endPos, lastStartPos, lastEndPos, startFloat, endFloat, lastStartFloat, lastEndFloat, min_length, min_markers, errorCount, lastErrorCount, lastWindowErrorCount, startPos2, endPos2, lastStartPos2, lastEndPos2, startFloat2, endFloat2, lastStartFloat2, lastEndFloat2, min_length2, min_markers2, errorCount2, lastErrorCount2, lastWindowErrorCount2, chrom, indiv1, indiv2, names, SNPID, geneMap, physMap, pFile, realIBD2, realIBD2New, noMerge);
	}

	printed=false;

	if((!merged) && (startFloat2!=-1)){//check the stored segment, print if valid
		segMarkerLength = endPos2-startPos2;
		segCMLength = endFloat2-startFloat2;
		if( (segCMLength > min_length2) && (segMarkerLength > min_markers2)){
			pFile.printf("%s\t%s\t%s\t%i\t%i\tIBD%i\t%f\t%f\t%f\t%i\t%i\t%f\n", PersonLoopData::_allIndivs[indiv1]->getId(), PersonLoopData::_allIndivs[indiv2]->getId(), chrom, physMap[startPos2], physMap[endPos2], 2, startFloat2, endFloat2, segCMLength, segMarkerLength, int(errorCount2), errorCount2 / float(segMarkerLength));
			printed=true;
		}

		if(printed){
			lastStartPos2=startPos2;
			lastEndPos2=endPos2;
			lastStartFloat2=startFloat2;
			lastEndFloat2=endFloat2;
			lastErrorCount2=errorCount2;
			startPos2=-100;
			endPos2=-100;
			startFloat2=-1;
			errorCount2=0;
			 //Deal with the IBD1 contiguous with the processed IBD2 segment.
			handleIBD1PostIBD2<IO_TYPE>(startPos, endPos, lastStartPos, lastEndPos, startFloat, endFloat, lastStartFloat, lastEndFloat, min_length, min_markers, errorCount, lastErrorCount, lastWindowErrorCount, startPos2, endPos2, lastStartPos2, lastEndPos2, startFloat2, endFloat2, lastStartFloat2, lastEndFloat2, min_length2, min_markers2, errorCount2, lastErrorCount2, lastWindowErrorCount2, chrom, indiv1, indiv2, names, SNPID, geneMap, physMap, pFile, realIBD2, realIBD2New, noMerge);
		}
	}
	return printed;




}




//Determine the corresponding block of 64 and bit within that block for genotype data.
void findIndexAndShift(uint64_t marker, uint64_t indiv, uint64_t newMarkerBlocks, uint64_t &newIndex, uint64_t &newShift){
	newIndex=indiv*(newMarkerBlocks)+(2*(marker/64));//Finds the 64-bit region in the transposed data where the value (first value) fits. Second goes in the next block
	newShift = marker%64;//Finds how deep into that 64-bit regions the new bits need to be placed.
}

//Determine the corresponding block of 64 and bit within that block for missing genotype data in a difference sized matrix.
void findIndexAndShiftMiss(uint64_t marker, uint64_t indiv, uint64_t newMarkerBlocks, uint64_t &newIndex, uint64_t &newShift){
	newIndex=indiv * (newMarkerBlocks / 2) + (marker / 64);
	newShift = marker % 64;
}

//Find block corresponding to an individual and set of markers.
inline uint64_t findIndex(uint64_t markerBlock, uint64_t indiv, uint64_t newMarkerBlocks){
	return indiv*(2*newMarkerBlocks)+(2*markerBlock);
}

//Find blockcorresponding to an individual and set of markers in the missing data.
inline uint64_t findMissIndex(uint64_t markerBlock, uint64_t indiv, uint64_t newMarkerBlocks){
	return indiv*newMarkerBlocks+markerBlock;
}



//populate the IBD boolean and error boolean int with true or false depending on if there is valid IBD1 in the current window and if it had any errors in it.
//Populate the errorCount value with the number of errors. Necessary for merge, error density analyses, and skipping IBD2 analysis when unneeded. 
void compareWindowsIBD1(uint64_t hom11, uint64_t hom12, uint64_t hom21, uint64_t hom22, bool &error, int &errorCount, bool &IBD){
	uint64_t fullCompare = 	(hom11 & hom22) | (hom12 & hom21);//The bitwise comparison for iBD1
	if(fullCompare==0){//If error free, return.
		IBD=true;
		error=false;
		errorCount=0;
		return;
	}

	fullCompare=fullCompare-lowestSetBit(fullCompare);
	if(fullCompare==0){//If one error, return and set errros to true.
		IBD=true;
		error=true;
		errorCount=1;
		return;
	}

	if((fullCompare-lowestSetBit(fullCompare))==0){//If there is more than one error, IBD is false and errors are irrelevant. Error count is
		IBD=false;//TODO
		error=false;//TODO
		errorCount=2;
		return;
	}
	IBD=false;
	error=false;
	errorCount=3;
	return;
}

//populate the IBD boolean and error boolean int with true or false depending on if there is valid IBD1 in the current window and if it had any errors in it.
//Populate the errorCount value with the number of errors. Necessary for merge and error density analyses.
void compareWindowsIBD2(uint64_t hom11, uint64_t hom12, uint64_t hom21, uint64_t hom22, uint64_t missing1, uint64_t missing2, bool &error, int &errorCount, bool &IBD){
	//uint64_t fullCompare = 	(~(missing1 | missing2)) & ((hom11 & ~hom21) | (hom12 & ~hom22));//Archaic and incorrect
	uint64_t fullCompare = (~(missing1 | missing2)) & ((hom11 ^ hom21) | (hom12 ^ hom22));//The bitwise comparison for IBD2.
	if(fullCompare==0){//If no error, return true.
		IBD=true;
		error=false;
		errorCount=0;
		return;
	}

	fullCompare = fullCompare - lowestSetBit(fullCompare);//If one error, return true
	if(fullCompare==0){
		IBD=true;
		error=true;
		errorCount=1;
		return;
	}

	if((fullCompare-lowestSetBit(fullCompare))==0){//If two errors, for IBD2 we still return true.
		IBD=true;//TODO
		error=true;//TODO
		errorCount=2;
		return;
	}
	IBD=false;
	error=false;
	errorCount=3;
	return;
}

//Transpose the input genotypes into a matrix of individuals as rows and blocks of markers as columns to speed up our analyses. Populates transposedData with homozygosity bit sets, and missingData with bit sets of where data is missing.
void interleaveAndConvertData(uint64_t *transposedData, uint64_t *missingData, uint64_t blocks, uint64_t markers, uint64_t bytes, uint64_t individuals, std::vector<int> &starts, std::vector<int> &ends, std::vector<float> &geneMap) {

	printf("Organizing genotype data for analysis... "); 
	fflush(stdout);
	uint64_t dataBlocks = ((individuals+3)/4);//Breakdown of the 8 bit block counts based on the number of individuals in the input. Not the blocks in our own internal data format.
	//Need to check later if redundant with blocks. These were different quantities before.

	uint64_t newMarkerBlocks = 2*(ceil((float)markers/64));//two for each set of 64 markers. This is the marker block format we use in our internal data storage.
	uint64_t mask1 = 1;

	
	std::fill(transposedData, transposedData + individuals * newMarkerBlocks, 0 );
	std::fill(missingData, missingData + individuals * (newMarkerBlocks / 2), 0);

/*	for (uint64_t b = 0; b < newMarkerBlocks; b++) {

		for (uint64_t i = 0; i < individuals; i++) {
			transposedData[i*newMarkerBlocks+b]=0;//Preemptively zeroing our matrix so our bit addition does not create problems.
		}
	}


	for (uint64_t b = 0; b < newMarkerBlocks/2; b++) {
		for (uint64_t i = 0; i < individuals; i++) {
			missingData[i*(newMarkerBlocks/2)+b]=0;//Preemptively zeroing our matrix so our bit addition does not create problems.
		}
	}*/

	uint8_t *data = new uint8_t[dataBlocks];

	for (uint64_t m = 0; m < markers; m++) {
		PersonIO<PersonLoopData>::readGenoRow(data, (int)(dataBlocks));//Reading in one row of the .bed file input to be transposed.
		uint64_t *data64 = (uint64_t *)data;//Reformat to allow blocks of 64 bits to be extracted and processed at once from an otherwise 8 bit format.

		uint8_t *data8 = data;
		for (uint64_t b = 0; b < blocks; b++) {

			uint64_t individuals1;
			uint64_t individuals2;
			if (b < individuals / 64) {
				individuals1 = data64[b * 2];//First set of 32 2-bit individuals to be combined into two blocks of 64 later.
				individuals2 = data64[b * 2 + 1];//Second set of 32 2-bit individuals to be combined into two blocks of 64 later. TODO May not be a speedup in new data ordering.
			}
			else {
				//Handle if the first set of individuals does not fill out a set of 64 completely and needs to be supplemented with 0s.
				if(individuals%64<32){
					individuals1 = 0;
					individuals2 = 0;
					int count = 0;
					for (uint64_t x = ((individuals-1) / 4); x>=b*16; x--) {
						individuals1 = individuals1 << 8;
						individuals1 += data8[x];
						count++;
					}
				}
				else{//Handle if the second set of individuals does not fill out a set of 64 completely and needs to be supplemented with 0s.
					individuals1 = data64[b * 2];
					individuals2 = 0;
					int count = 0;
					for (uint64_t x = ((individuals-1) / 4); x>=b*16+8; x--) {
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
				if(indiv>=individuals)
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
				findIndexAndShift(m, indiv, newMarkerBlocks, newIndex, newShift);
				uint64_t hom1Shift = hom1<<newShift;
				transposedData[newIndex]+=hom1Shift;

				uint64_t hom2Shift = hom2<<newShift;
				transposedData[newIndex+1]+=hom2Shift;



				findIndexAndShiftMiss(m, indiv, newMarkerBlocks, newIndex, newShift);
				uint64_t missShift = missBits<<newShift;
				missingData[newIndex]+=missShift;

			}
			for(uint64_t indivShift=32; indivShift<64; indivShift++){
				uint64_t indiv = indivShift+(b * 64);
				if(indiv>=individuals)
					break;
				uint64_t bit1 = (individuals2 >> (indivShift * 2)) & mask1;
				uint64_t bit2 = ((individuals2 >> (indivShift * 2+1)) & mask1);

				uint64_t hom1 = bit1 & bit2;
				uint64_t hom2 = ((~bit1) & (~bit2)) & mask1;

				uint64_t missBits = bit1 & ~bit2 & mask1;

				findIndexAndShift(m, indiv, newMarkerBlocks, newIndex, newShift);

				uint64_t hom1Shift = hom1<<newShift;
				transposedData[newIndex]+=hom1Shift;


				uint64_t hom2Shift = hom2<<newShift;
				transposedData[newIndex+1]+=hom2Shift;
				findIndexAndShiftMiss(m, indiv, newMarkerBlocks, newIndex, newShift);
				uint64_t missShift = missBits<<newShift;
				missingData[newIndex]+=missShift;

			}

		}
	}
	delete[] data;
	printf("done.\n");
}



//Loop over the individuals and windows and perform the segment analysis.
template<class IO_TYPE>
void segmentAnalysis(uint64_t *transposedData, uint64_t *missingData, int numIndivs, int indBlocks, int numMarkers, int markBytes, int markBlocks, std::vector<int> &starts, std::vector<int> &ends, std::vector<std::string> &names, std::vector<std::string> &SNPID, std::vector<float> &geneMap, std::vector<int> &physMap, std::string filename, std::string extension, bool noMerge, const char* chrom, int min_markers, float min_length, int min_markers2, float min_length2, float errorThreshold, float errorThreshold2, bool ibd2, int numThreads){	


	//Transposes the input matrix and handles the 8->64 bit format conversion..
	interleaveAndConvertData(transposedData, missingData, indBlocks, numMarkers, markBytes, numIndivs, starts, ends, geneMap);


	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(numThreads); //Forces input specified thread number with default 1. 

	printf("Beginning segment detection with %i thread(s)...",numThreads);
	fflush(stdout);

	//Begin parallelization
#pragma omp parallel
	{
		std::string threadname;
		FileOrGZ<IO_TYPE> pFile;
		threadname = filename +"."+std::to_string((1+omp_get_thread_num()))+extension;
		bool success = pFile.open(threadname.c_str(), "w");
		if(!success){
			printf("\nERROR: could not open output VCF file %s!\n", threadname.c_str());
			perror("open");
			exit(1);
		}

		//Grouping loop contents per parallelization, even though the individual threads start before this point.
		//
		//NOTE: -100 is used as a placeholder for marker positions of invalid start and endpoints to prevent merge from being viable with them, as merge cannot bridge gaps between a real position and something that far before 0 regardless of circumstance.
#pragma omp for schedule(dynamic, 10)
		for(uint64_t indiv1=0; indiv1<numIndivs; indiv1++){
			for(uint64_t indiv2=indiv1+1; indiv2<numIndivs; indiv2++){

				int startPos = -100;//NOTE! THIS ALLOWS ERRORS IN THE VERY FIRST WINDOW! TODO
				float startFloat = -1;//geneMap[startPos];//genetic start position of ongoing segment.
				float errorCount = 0;//error count of ongoing segment.
				int lastErrorFree = -100;//position of last error free window. Used to backtrack segment endings to the last error free window.
				float lastErrorFreeFloat = -1;//genetic position of the last error free window. NOTE: not part of the "last" group of variables referring to the previous stored segment.


				//Variables for storing data of the previous segment. Necessary for merging.			
				int lastStartPos = -100;//avoids marker merge issues
				int lastEndPos = -100;
				float lastEndFloat = -1;
				float lastStartFloat=-1;
				float lastErrorCount = 0;
				int lastWindowErrorCount = 0;//Part of merge


				int windowErrorCount = 0;//Part of merge.


				//Used to force print short IBD1 segments when attached to a real IBD2 segment.
				bool realIBD2 = false; //For IBD2-linked ongoing segments.
				bool realIBD2New = false; //Temporary holder for IBD2-linked status of new IBD1 segments created when an IBD2 segment split an ongoing IBD1 segment.


				//IBD2 equivalents of the previously described variables.
				int startPos2 = -100;//NOTE! THIS ALLOWS ERRORS IN THE VERY FIRST WINDOW! TODO
				float startFloat2 = -1;//geneMap[startPos];
				float errorCount2 = 0;
				int lastErrorFree2 = -100;
				float lastErrorFreeFloat2 = -1;

				int lastStartPos2 = -100;//avoids marker merge issues
				int lastEndPos2 = -100;
				float lastEndFloat2 = -1;
				float lastStartFloat2=-1;
				float lastErrorCount2 = 0;
				int lastWindowErrorCount2 = 0;
				int windowErrorCount2 = 0;







				for(uint64_t markerBlock=0; markerBlock<markBlocks; markerBlock++){

					//Booleans for the status of the IBD conclusions for the current window comparison.
					bool isIBD1;
					bool isIBD2;


					//Booleans for error status of the current window
					bool hasError;
					bool hasError2;


					//Locations of the two windows for comparison in the binary data.
					uint64_t index1 = findIndex(markerBlock, indiv1, markBlocks);
					uint64_t index2 = findIndex(markerBlock, indiv2, markBlocks);
					uint64_t missIndex1 = findMissIndex(markerBlock,indiv1,markBlocks);
					uint64_t missIndex2 = findMissIndex(markerBlock,indiv2,markBlocks);


					//Finds the two homozygous sets for each individual.
					uint64_t tDataI1_1 = transposedData[index1];//first individual, first homozygous set
					uint64_t tDataI1_2 =  transposedData[index1+1];//first individual, second homozygous set
					uint64_t tDataI2_1 = transposedData[index2];//second individual, first homozygous set
					uint64_t tDataI2_2 =  transposedData[index2+1];//second individual, second homozygous set.



					compareWindowsIBD1(tDataI1_1, tDataI1_2, tDataI2_1, tDataI2_2, hasError, windowErrorCount, isIBD1);

					//Don't bother checking IBD2 if there are too many errors in IBD1.
					if(ibd2 && windowErrorCount<3){
						compareWindowsIBD2(tDataI1_1, tDataI1_2, tDataI2_1, tDataI2_2, missingData[missIndex1], missingData[missIndex2], hasError2, windowErrorCount2, isIBD2);
					}
					else{
						isIBD2=false;
						windowErrorCount2=3;
						hasError2=false;
					}


					if(ibd2){//The check for if IBD2 is enabled
						if(isIBD2){//The check for if the window is IBD2

							int endPos2 = ends[markerBlock];//Set the endpos if the current window is the end of a segment.

							if(hasError2 && startPos2 != -100){
								errorCount2+=windowErrorCount2;
								//Error rate in ongoing segment check.
								if((errorCount2 / (endPos2 - startPos2)) > errorThreshold2){
									//Terminate ongoing IBD2 segment. Tracks if a segment was actually produced with validMerge.
									bool validMerge = endSegment(startPos2, lastErrorFree2, lastStartPos2, lastEndPos2, startFloat2, lastErrorFreeFloat2, lastStartFloat2, lastEndFloat2, min_length2, min_markers2, errorCount2, lastErrorCount2, lastWindowErrorCount2, 2, chrom, indiv1, indiv2, names, SNPID, geneMap, physMap, pFile, realIBD2, realIBD2New, noMerge);
									lastWindowErrorCount = windowErrorCount;//Tracks this for merge purposes. 
									if(validMerge){//If a segment really did get produced and printed, deals with the stored and ongoing iBD1 segments based on the position of the IBD2 segment.
										handleIBD1PostIBD2(startPos, lastErrorFree, lastStartPos, lastEndPos, startFloat, lastErrorFreeFloat, lastStartFloat, lastEndFloat, min_length, min_markers, errorCount, lastErrorCount, lastWindowErrorCount, startPos2, lastErrorFree2, lastStartPos2, lastEndPos2, startFloat2, lastErrorFreeFloat2, lastStartFloat2, lastEndFloat2, min_length2, min_markers2, errorCount2, lastErrorCount2, lastWindowErrorCount2, chrom, indiv1, indiv2, names, SNPID, geneMap, physMap, pFile, realIBD2, realIBD2New, noMerge);
									}
								}

							}
							else if(!hasError2){//If error free, no checks to error density need to be made here.
								if(startPos2 == -100){
									startPos2 = starts[markerBlock];
									startFloat2 = geneMap[starts[markerBlock]];
									errorCount2 = 0;
								}
								lastErrorFree2 = ends[markerBlock];
								lastErrorFreeFloat2 = geneMap[lastErrorFree2];
							}

						}
						else if(startPos2!=-100){//If no valid IBD2 was found, attempt to terminate the segment.
							bool validMerge = endSegment(startPos2, lastErrorFree2, lastStartPos2, lastEndPos2, startFloat2, lastErrorFreeFloat2, lastStartFloat2,lastEndFloat2, min_length2, min_markers2, errorCount2, lastErrorCount2, lastWindowErrorCount2, 2, chrom, indiv1, indiv2, names, SNPID, geneMap, physMap, pFile, realIBD2, realIBD2New, noMerge);
							if(validMerge){
								handleIBD1PostIBD2(startPos, lastErrorFree, lastStartPos, lastEndPos, startFloat, lastErrorFreeFloat, lastStartFloat, lastEndFloat, min_length, min_markers, errorCount, lastErrorCount, lastWindowErrorCount, startPos2, lastErrorFree2, lastStartPos2, lastEndPos2, startFloat2, lastErrorFreeFloat2, lastStartFloat2, lastEndFloat2, min_length2, min_markers2, errorCount2, lastErrorCount2, lastWindowErrorCount2, chrom, indiv1, indiv2, names, SNPID, geneMap, physMap, pFile, realIBD2, realIBD2New, noMerge);
							}
							lastWindowErrorCount2 = windowErrorCount2;
						}
					}




					if(isIBD1){

						int endPos = ends[markerBlock];
						if(hasError && startPos != -100){
							errorCount++;
							if((errorCount / (endPos - startPos)) > errorThreshold){
								//Terminates IBD2 segment if the IBD1 segment error density underneath it is too high.
								if(isIBD2){
									bool validMerge = endSegment(startPos2, lastErrorFree2, lastStartPos2, lastEndPos2, startFloat2, lastErrorFreeFloat2, lastStartFloat2, lastEndFloat2, min_length2, min_markers2, errorCount2, lastErrorCount2, lastWindowErrorCount2, 2, chrom, indiv1, indiv2, names, SNPID, geneMap, physMap, pFile, realIBD2, realIBD2New, noMerge);
									lastWindowErrorCount = windowErrorCount;
									if(validMerge){
										handleIBD1PostIBD2(startPos, lastErrorFree, lastStartPos, lastEndPos, startFloat, lastErrorFreeFloat, lastStartFloat, lastEndFloat, min_length, min_markers, errorCount, lastErrorCount, lastWindowErrorCount, startPos2, lastErrorFree2, lastStartPos2, lastEndPos2, startFloat2, lastErrorFreeFloat2, lastStartFloat2, lastEndFloat2, min_length2, min_markers2, errorCount2, lastErrorCount2, lastWindowErrorCount2, chrom, indiv1, indiv2, names, SNPID, geneMap, physMap, pFile, realIBD2, realIBD2New, noMerge);
									}
								}
								if(startPos!=-100 && endPos!=-100){
									endSegment(startPos, lastErrorFree, lastStartPos, lastEndPos, startFloat, lastErrorFreeFloat, lastStartFloat, lastEndFloat, min_length, min_markers, errorCount, lastErrorCount, lastWindowErrorCount, 1, chrom, indiv1, indiv2, names, SNPID, geneMap, physMap, pFile, realIBD2, realIBD2New, noMerge);
									lastWindowErrorCount = windowErrorCount;
								}
							}

						}
						else if(!hasError){
							if(startPos == -100){
								startPos = starts[markerBlock];
								startFloat = geneMap[starts[markerBlock]];
								errorCount = 0;
							}
							lastErrorFree = ends[markerBlock];
							lastErrorFreeFloat = geneMap[lastErrorFree];
						}

					}
					else if(startPos!=-100){
						endSegment(startPos, lastErrorFree, lastStartPos, lastEndPos, startFloat, lastErrorFreeFloat, lastStartFloat,lastEndFloat, min_length, min_markers, errorCount, lastErrorCount, lastWindowErrorCount, 1, chrom, indiv1, indiv2, names, SNPID, geneMap, physMap, pFile, realIBD2, realIBD2New, noMerge);
						lastWindowErrorCount = windowErrorCount;//TODO choose if this should be placed inside endSegment to reduce redundant code.
					}
				}
				//Handle the end of chromosomes, as there is no end of IBD or error density increase to trigger these segments to print otherwise.
				if(startFloat!=-1 || lastStartFloat!=-1 || startFloat2!=-1 || lastStartFloat2!=-1){//TODO put back forced
					if(ibd2){
						//Attempt to print IBD2 at the end of the chromosome if there is a valid segment.
						forceEndSegment2(startPos2, lastErrorFree2, lastStartPos2, lastEndPos2, startFloat2, lastErrorFreeFloat2, lastStartFloat2, lastEndFloat2, min_length2, min_markers2, errorCount2, lastErrorCount2, lastWindowErrorCount2, chrom, indiv1, indiv2, names, SNPID, geneMap, physMap, pFile, realIBD2, realIBD2New, startPos, lastErrorFree, lastStartPos, lastEndPos, startFloat, lastErrorFreeFloat, lastStartFloat, lastEndFloat, min_length, min_markers, errorCount, lastErrorCount, lastWindowErrorCount, noMerge);
					}
					//Attempt to print IBD1 at the end of the chromosome if there is a valid segment.
					forceEndSegment(startPos, lastErrorFree, lastStartPos, lastEndPos, startFloat, lastErrorFreeFloat, lastStartFloat, lastEndFloat, min_length, min_markers, errorCount, lastErrorCount, lastWindowErrorCount, 1, chrom, indiv1, indiv2, names, SNPID, geneMap, physMap, pFile, realIBD2, realIBD2New, noMerge);
				}
			}
		}
		pFile.close();
	}
	printf("done.\n");
}


//Consistent Print Statement for usage.
void printUsageAndExit(){

	printf("REQUIRED ARGUMENTS:\n");
	printf("\t EITHER:\n");
	printf("\t First Three Arguments [bed file] [bim file] [fam file]\t\tSpecifies the plink format files for the data by specific name.\n");
	printf("\t\t\t\t\t Must be first 3 arguments.\n");
	printf("\t OR:\n");
	printf("\t -b [prefix] or -bfile [prefix]\t\tSpecifies the prefix to be used with prefix.bed, prefix.bim, and prefix.fam for the plink format input.\n");
	printf("\t\t\t\t\t Does not need to be first argument\n.");
	printf("\n");
	printf("OPTIONS:\n");
	printf("\t-mL or -min_l <value>\t\t Specify minimum length for acceptible segments to output.\n");
	printf("\t\t\t\t\t Defaults to 7 centimorgans.\n");
	printf("\t-er or -errorRate <value> \t specify acceptible error rate in a segment before considering it false.\n");
	printf("\t\t\t\t\t Defaults to .004 errors per marker\n");
	printf("\t -f or file <filename>\t\t Specify output file.\n");
	printf("\t\t\t\t\t Defaults to ibis<thread number>.seg and will output a separate output file for each thread.\n");
	printf("\t -m or -merge\t\t\t allows internal merging of segments over gaps caused by potential errors.");
	printf("\n");
	printf("\t -2 or -ibd2\t\t\t enable ibd2 analyses");
	printf("\n");
	printf("\t -mt <value>\t\t\t set minimum number of markers required for acceptible segments to output.\n");
	printf("\t\t\t\t\t Defaults to 500 markers\n");
	printf("\n");
	printf("\t -threads <value>\t\t set the number of threads available to IBIS for parallel processing.");
	printf("\t\t\t\t\t Defaults to 1\n");
	printf("\n");
	printf("\t -gzip \t\t\t\t have the program output gzipped segment files\n");
	printf("OUTPUT FORMAT:\n");
	printf("sample1\tsample2\tchrom\tphys_start_pos\tphys_end_pos\tIBD_type\tgen_start_pos\tgen_end_pos\tseg_length\tmarker_count\tdebug_data\tdebug_data\n");
	printf("\t\n");
	exit(1);


}

int main(int argc, char **argv) {

	const char* VERSION_NUMBER = "1.05";
	const char* RELEASE_DATE = "November 10, 2019";
	printf("IBIS Segment Caller!  v%s    (Released %s)\n\n", VERSION_NUMBER, RELEASE_DATE);

	uint64_t numIndivs, numMarkers;//counts of input quantities.
	std::vector<std::string> names;
	std::vector<std::string> SNPID;
	std::vector<float> geneMap;//Genetic map of input
	std::vector<int> physMap;//Physical map of input
	float min_length = 7.0, min_length2 = 2.0;//cM minimum lengths
	bool noFile = true;//Check if no input file is given.
	bool ibd2 = false;//Check if IBD2 is requested.
	bool noMerge = true;//Default of not using merge.
	bool bFileNamesGiven = false;//Toggle for input as one argument or 3.
	bool distForce = false;//If true, stops the program from converting the input to cM.
	bool gzip = false;//If True, gzips the output.
	float errorDensityThreshold = 0.004, errorDensityThreshold2 = 0.008;//Maximum allowed error rates.
	float min_markers = 500, min_markers2 = 195;//Marker minimums for segments.
	const char* chrom = NULL;//Checks for specified chromosome. Optional, but required if the input contains multiple chromosomes.
	int numThreads;//Input threadnumber.
	numThreads=0;
	std::string filename, extension;//For building the output file.
	char bfileNameBed[100],bfileNameBim[100],bfileNameFam[100];//locations for input filenames.


	printf("Viewing arguments...\n");
	fflush(stdout);
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
		else if (arg == "-m" || arg == "-merge"){
			noMerge=false;
			printf("%s - running with merge enabled\n",arg.c_str());
		}
		else if (arg == "-mL" || arg == "-min_l") {
			min_length = atof(argv[i + 1]);
			min_length2 = min_length/2.0;
			printf("%s - running with minimum IBD1 length %f and minimum IBD2 length %f\n",arg.c_str(),min_length,min_length2);
		}
		else if (arg == "-mL2" || arg == "-min_l2") {
			min_length2 =  atof(argv[i + 1]);
			printf("%s - running with minimum IBD2 length %f\n",arg.c_str(),min_length2);
		}
		else if( arg =="-er" || arg == "-errorRate"){
			errorDensityThreshold=atof(argv[i+1]);
			errorDensityThreshold2=errorDensityThreshold*2;
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
			noFile = false;
			printf("%s - setting output file %s<thread number>.seg\n",arg.c_str(), filename.c_str());

		}
		else if(arg=="-mt"){
			min_markers=atoi(argv[i+1]);
			min_markers2=min_markers/2.0;
			printf("%s - running with minimum IBD1 marker threshold of %f and minimum IBD2 marker threshold of %f\n",arg.c_str(), min_markers,min_markers2);
		} else if(arg=="-mt2"){
			min_markers2=atoi(argv[i+1]);
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
		else if(arg.at(0)=='-'){
			printUsageAndExit();
		}

	}
	if(noFile){
		std::cout<<"No input file given. Defaulting filenames to \"ibis<thread number>.seg\"\"\n";
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
	uint64_t markBytes = ceil(((float)numMarkers * 2) / (8 * sizeof(uint8_t)));
	uint64_t indBlocks = ceil((float)(numIndivs) / 64); //blocks are based on the number of individuals
	uint64_t markBlocks = ceil((float)(numMarkers) / 64);


	for(uint64_t x = 0; x < numMarkers; x++){
		geneMap.push_back(Marker::getMarker(x)->getMapPos());
		physMap.push_back(Marker::getMarker(x)->getPhysPos());
		chrom = Marker::getMarker(x)->getChromName();
	}

	printf("Total SNPs: %i\n", int(geneMap.size()));

	//Convert to cM if the input genetic map seems to be too short.
	if(!distForce){
		float len = geneMap.back() - geneMap.front();
		if(len < 6.0){
			printf("Chromosome map shorter than 6 units of genetic distance.\n Morgan input detected - Converting to centimorgans. (Prevent this by running with -force argument)\n");
			for(uint64_t x = 0; x< geneMap.size(); x++){
				geneMap[x]=geneMap[x]*100.0;
			}
		}
		printf("Total SNPs: %i\n", int(geneMap.size())); 
	}

	uint64_t *transposedData = new uint64_t[2 * numIndivs*markBlocks];
	uint64_t *missingData = new uint64_t[numIndivs*markBlocks];//TODO THIS IS SILLY!




	std::vector<int> starts,ends;


	printf("Defining Windows... ");
	fflush(stdout);

	for(uint64_t x = 0; x < numMarkers; x += 64){
		starts.push_back(x);
		ends.push_back(min(numMarkers-1, x + 63));
	}


	printf("done.\n");

	
	if(gzip){
		extension = ".seg.gz";
		segmentAnalysis<gzFile>(transposedData, missingData, numIndivs, indBlocks, numMarkers, markBytes, markBlocks, starts, ends, names, SNPID, geneMap, physMap, filename, extension, noMerge, chrom, min_markers, min_length, min_markers2, min_length2, errorThreshold, errorThreshold2, ibd2, numThreads);	
	}
	else{
		extension = ".seg";
		segmentAnalysis<FILE *>(transposedData, missingData, numIndivs, indBlocks, numMarkers, markBytes, markBlocks, starts, ends, names, SNPID, geneMap, physMap, filename, extension, noMerge, chrom, min_markers, min_length, min_markers2, min_length2, errorThreshold, errorThreshold2, ibd2, numThreads);	
	}


}
