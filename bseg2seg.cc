#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <unordered_map>
#include <vector>
#include <tuple>

using namespace std;

void printUsage(char **argv);
void readFam(vector<char *> &indexToId, char *famFile, bool noFamId);

struct Segment {
  int inds[2];
  uint16_t chromIdx;
  int startPhys;
  int endPhys;
  uint8_t ibdType;
  float startGenet;
  float endGenet;
  int numMarkers;
  int errorCount;
};

int main(int argc, char **argv) {
  if (argc < 3) {
    printUsage(argv);
    exit(1);
  }

  vector<char *> indexToId;


  // Read in each segment file, printing its contents in human readable format

  bool famRead = false;
  for (int arg = 2; arg < argc; arg++) {
    char *segFile = argv[arg];
    FILE *in = fopen(segFile, "r");
    if (!in) {
      fprintf(stderr, "ERROR: couldn't open %s for reading: skipping\n",
	      segFile);
      perror("open");
      continue;
    }

//    fprintf(stderr, "Processing file %s... ", segFile);

    // first byte is fixed as character 129:
    if (fgetc(in) != 129) {
      fprintf(stderr, "ERROR: file %s does not appear to be a valid binary segment file\n",
	      segFile);
      fprintf(stderr, "       skipping\n");
      continue;
    }

    size_t ret;

    // next entry is an integer count of the byte length of the chromosome
    // string
    int chromStrLen;
    ret = fread(&chromStrLen, sizeof(int), 1, in);
    if (ret != 1) {
      fprintf(stderr, "ERROR: unable to read chromosome string length from bseg\n");
      exit(5);
    }

    // next is the chromosome string
    vector<char *> chroms;
    char *fullChromStr = new char[chromStrLen];
    ret = fread(fullChromStr, sizeof(char), chromStrLen, in);
    if ((int) ret != chromStrLen) {
      fprintf(stderr, "ERROR: unable to read chromosome string from bseg\n");
      exit(5);
    }

    // now break apart the chromosome string into each chromosome
    chroms.push_back(fullChromStr);
    for(int i = 1; i < chromStrLen - 1; i++) { // -1: last char should be '\0'
      if (fullChromStr[i] == '\0')
	chroms.push_back( fullChromStr + i + 1 );
    }

    // next: a boolean (stored in a byte) indicating whether to omit family ids
    uint8_t noFamId;
    ret = fread(&noFamId, sizeof(uint8_t), 1, in);
    if (ret != 1) {
      fprintf(stderr, "ERROR: unable to read the no family id field from bseg\n");
      exit(5);
    }

    if (!famRead) {
      readFam(indexToId, /*famFile=*/ argv[1], noFamId);
      famRead = true;
    }

    Segment seg;
    while (fread(&seg.inds[0], sizeof(int), 1, in) != 0) {
      size_t allRet = fread(&seg.inds[1], sizeof(int), 1, in);

      ret = fread(&seg.chromIdx, sizeof(uint16_t), 1, in);
      allRet = allRet && ret;

      ret = fread(&seg.startPhys, sizeof(int), 1, in);
      allRet = allRet && ret;

      ret = fread(&seg.endPhys, sizeof(int), 1, in);
      allRet = allRet && ret;

      ret = fread(&seg.ibdType, sizeof(uint8_t), 1, in);
      allRet = allRet && ret;

      ret = fread(&seg.startGenet, sizeof(float), 1, in);
      allRet = allRet && ret;

      ret = fread(&seg.endGenet, sizeof(float), 1, in);
      allRet = allRet && ret;

      ret = fread(&seg.numMarkers, sizeof(int), 1, in);
      allRet = allRet && ret;

      ret = fread(&seg.errorCount, sizeof(int), 1, in);
      allRet = allRet && ret;

      if (allRet != 1) {
	fprintf(stderr, "ERROR: unable to read entry in middle of bseg file: was ibis interrupted?\n");
	exit(5);
      }

      printf("%s\t%s\t%s\t%d\t%d\tIBD%d\t%f\t%f\t%f\t%d\t%d\t%f\n",
	     indexToId[ seg.inds[0] ], indexToId[ seg.inds[1] ],
	     chroms[ seg.chromIdx ], seg.startPhys, seg.endPhys, seg.ibdType,
	     seg.startGenet, seg.endGenet, (seg.endGenet - seg.startGenet),
	     seg.numMarkers, seg.errorCount, (float) seg.errorCount / seg.numMarkers);
    }

//    fprintf(stderr, "done.\n");
  }
}

void printUsage(char **argv) {
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "%s [fam file] [bseg files ...]\n",
	  argv[0]);
}

void readFam(vector<char *> &indexToId, char *famFile, bool noFamId) {
  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  if (buffer == NULL) {
    fprintf(stderr, "ERROR: out of memory");
    exit(5);
  }

  FILE *in = fopen(famFile, "r");
  if (!in) {
    fprintf(stderr, "ERROR: couldn't open fam file %s for reading\n", famFile);
    fprintf(stderr, "       required for processing, so won't proceed\n");
    perror("open");
    exit(5);
  }

  int curSampIdx = 0;
  while (getline(&buffer, &bytesRead, in) >= 0) {
    const char *delim = " \t\n";
    char *saveptr;

    char *famId = strtok_r(buffer, delim, &saveptr);
    char *sampleId = strtok_r(NULL, delim, &saveptr);

    // confirm there are 4 more fields
    for(int i = 0; i < 4; i++) {
      char *dontcare = strtok_r(NULL, delim, &saveptr);
      if (dontcare == NULL) {
	fprintf(stderr, "ERROR: line %d of fam file contains fewer than 6 fields\n",
		curSampIdx + 1);
	fprintf(stderr, "       won't proceed\n");
	exit(1);
      }
    }

    char *dontcare = strtok_r(NULL, delim, &saveptr);
    if (dontcare != NULL) {
      fprintf(stderr, "ERROR: line %d of fam file contains more than 6 fields\n",
	      curSampIdx + 1);
      fprintf(stderr, "       won't proceed\n");
      exit(1);
    }

    // good entry: store
    char *id;
    if (noFamId) {
      id = new char[ strlen(sampleId) + 1 ]; // +1 for '\0'
      sprintf(id, "%s", sampleId);
    }
    else {
      id = new char[ strlen(famId) + strlen(sampleId) + 2 ]; // +2 ':', '\0'
      sprintf(id, "%s:%s", famId, sampleId);
    }
    indexToId.push_back( id );
    curSampIdx++;
  }

  fclose(in);

  free(buffer);
}
