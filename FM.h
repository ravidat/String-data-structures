/*  FM-Index - Text Index
 *  Copyright (C) 2011  Matthias Petri
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>. */
 
#ifndef FM_H
#define	FM_H

#include <stdarg.h>

#include "util.h"

/* libcds includes */
#include "libcdsBasics.h"
#include "BitString.h"
#include "BitSequence.h"
#include "BitSequenceRG.h"
#include "Sequence.h"
#include "WaveletTreeNoptrs.h"
#include "Mapper.h"
#include "BitSequenceBuilder.h"

using namespace std;
using namespace cds_utils;
using namespace cds_static;

#define DEFAULT_SAMPLERATE      64
#define RRR_SAMPLERATE			20

/*Struct for a single interval [sp,ep]*/
struct SR {
  uint32_t sp;
  uint32_t ep;
  SR(): sp(0), ep(0){};
  SR(uint32_t _sp, uint32_t _ep): sp(_sp), ep(_ep){}; /*Constructor*/
  void setRange(uint32_t _sp, uint32_t _ep){this->sp = _sp; this->ep = _ep;}; /*Set a new range*/
  uint32_t getSp() {return sp;}; /*Returns left interval end*/
  uint32_t getEp() {return ep;}; /*Returns right interval end*/
  void setEp(uint32_t _ep) {this->ep = _ep;}; /*Set a new ep*/
  void setSp(uint32_t _sp) {this->sp = _sp;}; /*Set a new sp*/ 
  uint32_t getSize() {return this->ep - this->sp + 1;}; /*Returns number of elements in range*/
  virtual ~SR(){};
}; 

/*Struct for a double interval [sp,ep] and [spR,epR]*/
struct SR2 {
  SR sr;
  SR srR;
  uint8_t valid;
  SR2():sr(),srR(),valid(0){}; /*Empty constructor*/
  SR2(SR _sr, SR _srR):sr(_sr),srR(_srR),valid(0){}; /*Constructor*/
  SR2(uint32_t sp, uint32_t ep, uint32_t spR, uint32_t epR):sr(sp,ep),srR(spR,epR),valid(0){}; /*Constructor that creates both SR instances from 2 given ranges*/
  SR getSR(){return this->sr;}; /*Returns the suffix range for the regular string*/
  SR getSRr(){return this->srR;}; /*Returns the suffix range for the reverse string*/
  void setSR(SR _sr){this->sr=_sr;}; /*Sets the suffix range for the regular string*/
  void setSRr(SR _srR){this->srR=_srR;}; /*Sets the suffix range for the reverse string*/
  uint32_t getSize(){return this->sr.getSize();}; /*Returns the size of the range of the first element (the second has the same size)*/
  bool isLegalRange(){return ((sr.getSp() <= sr.getEp()) && (srR.getSp() <= srR.getEp()));}; /*Returns TRUE iff the range of SR and SR reverse is legal*/
  virtual ~SR2(){};
};

class FM {
public:
    FM(uint8_t* T,uint32_t n,uint32_t samplerate);
    void build(uint8_t* T,uint32_t n,uint32_t samplerate);
    static FM* load(char* filename);
    int32_t save(char* filename);
    uint8_t* remap0(uint8_t* T,uint32_t n);
    uint32_t count(uint8_t* pattern,uint32_t m);
    uint32_t* locate(uint8_t* pattern,uint32_t m,uint32_t* matches);
    uint32_t getSize();
    uint8_t* extract(uint32_t start,uint32_t stop);
    uint8_t* reconstructText(uint32_t* n);
    float getSizeN();
    virtual ~FM();
    /******ADDED FUNCTIONS******/
    SR2 reverseSR_SRr(SR2 sr2);
    size_t smallerCharSuffixes(uint8_t c, uint32_t bound, WaveletTreeNoptrs *T);
    SR2 StartBackwardSearch(uint8_t* pattern, SR2 sr2, int32_t sp, int32_t ep); 
    SR2 getNewSR2(SR2 sr2, uint8_t c, WaveletTreeNoptrs *T); 
    SR getNewBackSR(SR sr, uint8_t c, WaveletTreeNoptrs *T);
    SR2 StartForwardSearch(uint8_t* pattern, SR2 sr2, int32_t sp, int32_t ep); 
    //SR2* forwardSearch(uint8_t* X, uint8_t* pattern, SR2 sr2, uint8_t c, uint32_t sp, uint32_t ep);
    uint8_t* reverse(uint8_t* X, uint32_t n);
    std::vector<SR2> locate1err(uint8_t* pattern,uint32_t m);
    std::vector<SR2> locate2err(uint8_t* pattern,uint32_t m);
    uint32_t* getIndexes(SR sr, uint32_t *matches);
    size_t getRank(WaveletTreeNoptrs *T, uint8_t c, int32_t range);
    void print(SR2 sr2);
    void printPtr(SR2 *sr2);


public:
	static void info(const char *format,...)
	{
		if(FM::verbose == 1) {
			va_list vargs;
			va_start (vargs, format);
			vfprintf (stderr, format, vargs);
			fprintf (stderr, "\n");
		}
	}
	static int verbose;
private:
    FM();
private:
  uint32_t sigma;
  uint32_t samplerate;
  int32_t I;
  uint32_t n;
  uint32_t C[size_uchar+1]; /*C[i] is equal to the number of characters in T smaller than i */
  uint8_t remap[size_uchar];
  uint8_t* remap_reverse;
  uint32_t* suffixes;
  uint32_t* positions;
  BitSequence* sampled;
  WaveletTreeNoptrs *T_bwt;
  WaveletTreeNoptrs *Tr_bwt; /*Added reverse tree*//*TODO: free memory*/
  Sequence *N;
  Sequence *M;
};

#endif	/* FM_H */

