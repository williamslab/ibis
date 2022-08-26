#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <unordered_map>
#include <vector>
#include <tuple>

using namespace std;

struct HashString {
  size_t operator() (const char *s) const {
    // got the hash from stack overflow and changed macros to constants
    // https://stackoverflow.com/questions/8317508/hash-function-for-a-string
    const int A = 54059; /* a prime */
    const int B = 76963; /* another prime */
    //const int C = 86969; /* yet another prime */
    const int FIRSTH = 37; /* also prime */

    size_t h = FIRSTH;
    while (*s) {
      h = (h * A) ^ (s[0] * B);
      s++;
    }

    return h; // or return h % C;
  }
};

struct EqualString {
  bool operator() (const char *lhs, const char *rhs) const {
    return !strcmp(lhs, rhs);
  }
};

struct SampIdIdx {
  unordered_map<char*, int, HashString, EqualString> idTo;
  vector<char *> toId;
};

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


void printUsage(char **argv);
void readFam(SampIdIdx &idxs, char *famFile, bool noFamId);
uint64_t readSeg(FILE *in, char *segFile, char *buffer, size_t &bytesRead,
		 bool &famRead, char *famFile, SampIdIdx &idxs,
		 unordered_map<int, tuple<int, float, float>> *&pairIBDsums);
uint64_t readBseg(FILE *in, char *segFile, bool &famRead, char *famFile,
		  SampIdIdx &idxs,
		  unordered_map<int, tuple<int, float, float>> *&pairIBDsums);
void addSeg(unordered_map<int, tuple<int, float, float>> *&pairIBDsums,
	    int sampIdx[2], uint8_t ibdType, float genetLength);


int main(int argc, char **argv) {
  if (argc < 4) {
    printUsage(argv);
    exit(1);
  }

  char *endptr;
  char *mapLengthStr = argv[1];
  float mapLength = strtof(mapLengthStr, &endptr);
  if (errno != 0 || *endptr != '\0') {
    fprintf(stderr, "ERROR: unable to convert %s to floating point\n",
	    mapLengthStr);
    if (errno != 0)
      perror("strtof");

    fprintf(stderr, "\n");
    printUsage(argv);
    exit(6);
  }

  // sample ids to index map
  SampIdIdx idxs;

  // Logically a 2d "array", with both indexes being integer indexes s1 and s2
  // (from idToIndex)
  // We require that s1 < s2
  // The value stores s1 and s2s:
  //   (0) number of IBD segments
  //   (1) the cM length of their IBD1 segments
  //   (2) the cM length of their IBD2 segments
  unordered_map< int, tuple<int, float, float>> *pairIBDsums = NULL;

  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  if (buffer == NULL) {
    fprintf(stderr, "ERROR: out of memory");
    exit(5);
  }


  // read in each segment file, storing the sums of the segment lengths

  bool famRead = false;
  for (int arg = 3; arg < argc; arg++) {
    char *segFile = argv[arg];
    FILE *in = fopen(segFile, "r");
    if (!in) {
      fprintf(stderr, "ERROR: couldn't open %s for reading: skipping\n",
	      segFile);
      perror("open");
      continue;
    }

    fprintf(stderr, "Processing file %s... ", segFile);

    uint64_t nSegs;

    // first byte of bseg file is 129 (won't appear in plain text seg)
    int c = fgetc(in);
    if (c == 129) {
      nSegs = readBseg(in, segFile, famRead, /*famFile=*/ argv[2], idxs,
		       pairIBDsums);
    }
    else {
      ungetc(c, in); // put first byte back

      nSegs = readSeg(in, segFile, buffer, bytesRead, famRead,
		      /*famFile=*/ argv[2], idxs, pairIBDsums);
    }
    fprintf(stderr, "\rProcessing file %s... %ld segments total... done.  \n",
	    segFile, nSegs);

  }

  fprintf(stderr, "Done processing segments.\n\n");

  // LAST: print coef file

  // degree limits
  float degThresh[8];
  degThresh[0] = 0.475;
  degThresh[1] = 0.176776695;
  degThresh[2] = 0.088388348;
  degThresh[3] = 0.044194174;
  degThresh[4] = 0.022097087;
  degThresh[5] = 0.011048543;
  degThresh[6] = 0.005524272;
  degThresh[7] = 0.002762136;

  fprintf(stderr, "Printing... ");
  // header
  printf("Individual1\tIndividual2\tKinship_Coefficient\tIBD2_Fraction\tSegment_Count\tDegree\n");

  int numSamples = idxs.idTo.size();
  for(int s0 = 0; s0 < numSamples; s0++) {
    char *s0Id = idxs.toId[ s0 ];
    for(auto it = pairIBDsums[s0].begin(); it != pairIBDsums[s0].end(); it++) {
      char *s1Id = idxs.toId[ it->first ];
      auto &theTuple = it->second;

      int segCount = get<0>(theTuple);
      float ibd1Length = get<1>(theTuple);
      float ibd2Length = get<2>(theTuple);
      float phi = (.25 * ibd1Length + .5 * ibd2Length) / mapLength;
      // TODO: option to set additional factor
      float addPhi = 0.001380f;
      phi += addPhi;
      int degree;
      if (phi < degThresh[7])
	degree = -1;
      else {
	for(degree = 0; degree < 8; degree++) {
	  if (phi > degThresh[degree])
	    break;
	}
      }
      printf("%s\t%s\t%.6f\t%.6f\t%d\t%d\n",
	     s0Id, s1Id, phi, ibd2Length / mapLength, segCount, degree);
    }
  }

  fprintf(stderr, "done.\n");
}

void printUsage(char **argv) {
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "%s [total map length cM] [fam file] [seg/bseg files ...]\n",
	  argv[0]);
}

uint64_t readSeg(FILE *in, char *segFile, char *buffer, size_t &bytesRead,
		 bool &famRead, char *famFile, SampIdIdx &idxs,
		 unordered_map<int, tuple<int, float, float>> *&pairIBDsums) {
  char *endptr;

  uint64_t line = 1;
  while (getline(&buffer, &bytesRead, in) >= 0) {
    if (line % 1000000 == 0) {
      fprintf(stderr, "\rProcessing file %s... at line %ld million",
	      segFile, line / 1000000);
    }
    const char *delim = " \t\n";
    char *saveptr;
    char *samps[2];
    int sampIdx[2];

    samps[0] = strtok_r(buffer, delim, &saveptr);
    samps[1] = strtok_r(NULL, delim, &saveptr);

    if (!famRead) {
      bool noFamId = true; // assume we won't be concatenating the family id
      for(int c = 0; samps[0][c] != '\0'; c++)
	if (samps[0][c] == ':')
	  noFamId = false; // have ':' family delimiter: want family id
      readFam(idxs, famFile, noFamId);
      famRead = true;

      // now make space for the sums:
      int numSamples = idxs.idTo.size();
      // TODO: resize?
      pairIBDsums =
		  new unordered_map< int, tuple<int, float, float>>[numSamples];
    }

    for(int i = 0; i < 2; i++) {
      auto entry = idxs.idTo.find(samps[i]);
      if (entry == idxs.idTo.end()) {
	fprintf(stderr, "ERROR: id %s in seg file %s not in fam file (line %ld of seg)\n",
		samps[i], segFile, line);
	fprintf(stderr, "       won't proceed\n");
	exit(3);
      }
      sampIdx[i] = entry->second;
    }

    strtok_r(NULL, delim, &saveptr); /*chrom*/
    strtok_r(NULL, delim, &saveptr); /*physStart*/
    strtok_r(NULL, delim, &saveptr); /*physEnd*/
    char *ibdTypeStr = strtok_r(NULL, delim, &saveptr);
    strtok_r(NULL, delim, &saveptr); /*genetStart*/
    strtok_r(NULL, delim, &saveptr); /*genetEnd*/
    char *genetLengthStr = strtok_r(NULL, delim, &saveptr);

    // process IBD type
    uint8_t ibdType = 0;
    if (strcmp(ibdTypeStr, "IBD1") == 0)
      ibdType = 1;
    else if (strcmp(ibdTypeStr, "IBD2") == 0)
      ibdType = 2;
    else {
      fprintf(stderr, "ERROR: IBD type of segment is %s, which the code can't process\n",
	      ibdTypeStr);
      fprintf(stderr, "       (line %ld of %s)", line, segFile);
      fprintf(stderr, "       won't proceed\n");
      exit(2);
    }

    // get floating point length
    float genetLength = strtof(genetLengthStr, &endptr);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: unable to convert genetic length %s to floating point\n",
	      genetLengthStr);
      fprintf(stderr, "       (line %ld of %s)", line, segFile);
      fprintf(stderr, "       won't proceed\n");
      if (errno != 0)
	perror("strtof");
      exit(6);
    }

    assert(sampIdx[0] < sampIdx[1]);

    addSeg(pairIBDsums, sampIdx, ibdType, genetLength);

    line++;
  }

  return line - 1;
}

uint64_t readBseg(FILE *in, char *segFile, bool &famRead, char *famFile,
		  SampIdIdx &idxs,
		  unordered_map<int, tuple<int, float, float>> *&pairIBDsums) {
  // already read the first byte

  // next entry is an integer count of the byte length of the chromosome
  // string
  int chromStrLen;
  size_t ret = fread(&chromStrLen, sizeof(int), 1, in);
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

  // no need to break apart the chromosome string into each chromosome: we
  // don't use these strings for calculating the coefficients 

  // next: a boolean (stored in a byte) indicating whether to omit family ids
  uint8_t noFamId;
  ret = fread(&noFamId, sizeof(uint8_t), 1, in);
  if (ret != 1) {
    fprintf(stderr, "ERROR: unable to read the no family id field from bseg\n");
    exit(5);
  }

  if (!famRead) {
    readFam(idxs, famFile, noFamId);
    famRead = true;

    // now make space for the sums:
    int numSamples = idxs.idTo.size();
    // TODO: resize?
    pairIBDsums = new unordered_map< int, tuple<int, float, float>>[numSamples];
  }

  uint64_t nSegs = 0;
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

    // length of segment:
    float genetLength = seg.endGenet - seg.startGenet;

    assert(seg.inds[0] < seg.inds[1]);

    addSeg(pairIBDsums, seg.inds, seg.ibdType, genetLength);

    nSegs++;
  }

  return nSegs;
}

void readFam(SampIdIdx &idxs, char *famFile, bool noFamId) {
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

    // good entry: add to map:
    char *id;
    if (noFamId) {
      id = new char[ strlen(sampleId) + 1 ]; // +1 for '\0'
      sprintf(id, "%s", sampleId);
    }
    else {
      id = new char[ strlen(famId) + strlen(sampleId) + 2 ]; // +2 ':', '\0'
      sprintf(id, "%s:%s", famId, sampleId);
    }

    if (idxs.idTo.find(id) != idxs.idTo.end()) {
      fprintf(stderr, "ERROR: fam file contains multiple entries of id %s (line %d)\n",
	      id, curSampIdx+1);
      fprintf(stderr, "       won't proceed\n");
      exit(2);
    }
    idxs.idTo[id] = curSampIdx;
    idxs.toId.push_back( id );
    curSampIdx++;
  }

  fclose(in);

  free(buffer);
}

void addSeg(unordered_map<int, tuple<int, float, float>> *&pairIBDsums,
	    int sampIdx[2], uint8_t ibdType, float genetLength) {
  auto mapEntry = pairIBDsums[ sampIdx[0] ].find( sampIdx[1] );
  if (mapEntry == pairIBDsums[ sampIdx[0] ].end()) {
    tuple<int,float,float> theTuple = make_tuple(1, 0.0f, 0.0f);
    switch (ibdType) {
      case 1:
	get<1>(theTuple) = genetLength;
	break;
      case 2:
	get<2>(theTuple) = genetLength;
	break;
    }

    pairIBDsums[ sampIdx[0] ][ sampIdx[1] ] = theTuple;
  }
  else {
    get<0>( mapEntry->second )++;
    switch (ibdType) {
      case 1:
	get<1>( mapEntry->second ) += genetLength;
	break;
      case 2:
	get<2>( mapEntry->second ) += genetLength;
	break;
    }
  }
}
