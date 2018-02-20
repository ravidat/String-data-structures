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

#include "FM.h"
#include "util.h"
#include "divsufsort.h"
#include <math.h> /*Added*/
#include <algorithm>
 int FM::verbose = 0;


 FM::FM(uint8_t* T,uint32_t N,uint32_t samplerate = DEFAULT_SAMPLERATE) {
    this->samplerate = samplerate;
    this->n = N;
    
    /* 0 terminate */
    if(T[N-1] != 0) {
        T = (uint8_t*) safe_realloc(T,(N+1) * sizeof(uint8_t));
        T[N] = 0;
        this->n++;
    } 
    build(T,n,samplerate);
}

FM::FM() {
    this->n = 0;
    this->sigma = 0;
    this->I = 0;
    this->remap_reverse = NULL;
    this->T_bwt = NULL;
    this->Tr_bwt = NULL; /*Added initialize*/
    this->sampled = NULL;
    this->suffixes = NULL;
}

FM::~FM() {
    free(remap_reverse);
    free(suffixes);
    free(positions);
    delete T_bwt;
    delete Tr_bwt;
    delete sampled;
}

uint32_t
FM::getSize() {
    uint32_t bytes = 0;
    
    bytes += sizeof(this->n);
    bytes += sizeof(this->samplerate);
    bytes += sizeof(this->sigma);
    bytes += sizeof(this->I);
    bytes += sizeof(this->remap);
    bytes += sizeof(this->C);
    bytes += this->sigma * sizeof(uint8_t); /* remap_reverse */
    bytes += ((n/samplerate)+1) * sizeof(uint32_t); /* suffixes */
	bytes += ((n/samplerate)+2) * sizeof(uint32_t); /* positions */
    bytes += this->sampled->getSize();
    bytes += this->T_bwt->getSize();
    
    return bytes;
}

float
FM::getSizeN() {
    uint32_t bytes = getSize();
    return (float)(bytes)/(float)(n);
}

uint8_t*
FM::remap0(uint8_t* T,uint32_t n) {
    uint8_t* X;
    uint32_t i,j,size=0;
    uint32_t freqs[size_uchar]; /*size_uchar = number of possible characters*/
    
    for(i=0;i<size_uchar;i++) freqs[i]=0; /*Zero all elements in freqs array*/
    for(i=0;i<n;i++) if(freqs[T[i]]++==0) size++; /*At the end of the loop, size = Sigma size of T (The whole string)*/
    
    this->sigma=size;
    
    // remap alphabet
    if (freqs[0]>1) {i=1;sigma++;} //test if some character of T is zero, we already know that text[n-1]='\0'
    else i=0;

    remap_reverse = (uint8_t*) malloc(size*sizeof(uint8_t));
    for(j=0;j<size_uchar;j++) {
      if(freqs[j]!=0) { /*If the char with the value of j appears in T*/
        remap[j]=i;
        remap_reverse[i++]=j;
    }
}
    // remap text
X = (uint8_t*) malloc(n * sizeof(uint8_t));
    for(i=0;i<n-1;i++) // the last character must be zero
        X[i]=remap[T[i]];  
    return X;
}

/*Creating (for each of these, also the reverse is created): the X string, bwt array, suffix array, Sampled suffix array, wavelet tree*/
void
FM::build(uint8_t* T,uint32_t n,uint32_t samplerate) {
    uint8_t *X,/*reverse*/ *Xr;
    uint8_t *X_bwt, /*reverse*/ *Xr_bwt;
    int32_t *SA, /*reverse*/ *SAr;
    uint32_t i,prev,tmp,start,stop;
    float elapsed;
    
    start = gettime();

    info("building index.");
    
    /* remap if 0 in text */
    info("- remapping alphabet.");
    X = remap0(T,n);
    Xr = reverse(X,n); /*TODO: implement function and destructor*/
    free(T);
    
    /*C[i] is equal to the number of characters in T smaller than i */
    /* create cumulative counts */
    info("- creating cumulative counts C[].");
    for (i=0;i<size_uchar+1;i++) {C[i]=0;} /*Initialize all C elements as zero*/
    for (i=0;i<n;++i) C[X[i]]++;
        prev=C[0];C[0]=0;

    for (i=1;i<size_uchar+1;i++) {
      tmp = C[i];
      C[i]=C[i-1]+prev;
      prev = tmp;
  }
    /* perform k-BWT */
  info("- performing bwt.");
  SA = (int32_t*) safe_malloc( n * sizeof(int32_t)  );
    if( divsufsort(X,SA,n) != 0 ) { /*Sort the suffix array SA, without memory compression*/
  fatal("error divsufsort");
}

    /*perform k-BWT for the reverse T*/
info("- performing bwt reverse.");
SAr = (int32_t*) safe_malloc( n * sizeof(int32_t)  );
    if( divsufsort(Xr,SAr,n) != 0 ) { /*Sort the suffix array SA reverse, without memory compression*/
fatal("error divsufsort reverse");
}

    /*Creating the SSA (sampled suffix array) 
    We should not create an SSA for the reverse string, as it is too memory consuming*/
    /* sample SA for locate() */
info("- sample SA locations.");
    suffixes = (uint32_t*) safe_malloc( ((n/samplerate)+1) * sizeof(uint32_t)); /*suffixes will hold the SSA (containing n/samplerate strings)*/
    BitString B(n); /*Will hold 1 if i%samplerate = 1, 0 otherwise*/
tmp = 0;
for(i=0;i<n;i++) {
    if( SA[i] % samplerate == 0) {
        suffixes[tmp] = SA[i];
        B.setBit(i,true);
        tmp++;
    } else B.setBit(i,false);
}
    /* enable rank on context vector */
this->sampled = new BitSequenceRRR(B,RRR_SAMPLERATE);

	/* sample SA for display() */
positions = (uint32_t*) safe_malloc( ((n/samplerate)+2) * sizeof(uint32_t));
for (i=0;i<this->n;i++) {
    if (SA[i] % samplerate == 0) this->positions[SA[i]/samplerate] = i;
}
positions[(this->n-1)/samplerate+1] = positions[0];

info("- creating bwt output.");
X_bwt = (uint8_t*) safe_malloc( n * sizeof(uint8_t)  );
for(i=0;i<n;i++) {
    if(SA[i]==0) { 
        X_bwt[i] = X[n-1];
        this->I = i;
    } else X_bwt[i] = X[SA[i]-1];
}
free(SA);
info("- creating bwt reverse output.");
Xr_bwt = (uint8_t*) safe_malloc( n * sizeof(uint8_t)  );

for(i=0;i<n;i++) {
    if(SAr[i]==0) { 
        Xr_bwt[i] = Xr[n-1];
    } 
    else {
        Xr_bwt[i] = Xr[SAr[i]-1];
    }
}
    free(SAr);
    free(Xr);

info("- create RRR wavelet tree over bwt.");
MapperNone * map = new MapperNone();
BitSequenceBuilder * bsb = new BitSequenceBuilderRRR(RRR_SAMPLERATE);
T_bwt = new WaveletTreeNoptrs((uint32_t*)X_bwt,n,sizeof(uint8_t)*8,bsb,map,true);

info("- create RRR wavelet tree reverse over bwt.");
MapperNone * mapr = new MapperNone();
BitSequenceBuilder * bsbr = new BitSequenceBuilderRRR(RRR_SAMPLERATE);
Tr_bwt = new WaveletTreeNoptrs((uint32_t*)Xr_bwt,n,sizeof(uint8_t)*8,bsbr,mapr,true);

stop = gettime();
elapsed = (float)(stop-start)/1000000;

    /* build aux data */
info("build FM-Index done. (%.3f sec)",elapsed);

uint32_t bytes;
info("space usage:");
bytes = sigma * sizeof(uint8_t);
info("- remap_reverse: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
bytes = sizeof(this->C);
info("- C: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
bytes = ((n/samplerate)+1) * sizeof(uint32_t);
info("- Suffixes: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
bytes = ((n/samplerate)+2) * sizeof(uint32_t);
info("- Positions: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
bytes = sampled->getSize();
info("- Sampled: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
bytes = T_bwt->getSize();
info("- T_bwt: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100);
bytes = Tr_bwt->getSize();     /*Added for the Xr*/
info("- Tr_bwt: %d bytes (%.2f\%)",bytes,(float)bytes/getSize()*100); 
info("input Size n = %lu bytes\n",this->n);
info("index Size = %lu bytes (%.2f n)",getSize(),getSizeN());
}

/*Save */
int32_t
FM::save(char* filename) {
    std::ofstream f;
    f.open(filename,std::ios::out | std::ios::binary); 
    
    info("writing FM Index to file '%s'",filename);
    if(f.is_open()) {
        f.write(reinterpret_cast<char*>(&samplerate),sizeof(uint32_t));
        f.write(reinterpret_cast<char*>(&sigma),sizeof(uint32_t));
        f.write(reinterpret_cast<char*>(&I),sizeof(uint32_t));
        f.write(reinterpret_cast<char*>(&n),sizeof(uint32_t));
        f.write(reinterpret_cast<char*>(C),sizeof(uint32_t)*(size_uchar+1));
        f.write(reinterpret_cast<char*>(remap),sizeof(uint8_t)*size_uchar);
        f.write(reinterpret_cast<char*>(remap_reverse),sizeof(uint8_t)*sigma);
        f.write(reinterpret_cast<char*>(suffixes),sizeof(uint32_t)*((n/samplerate)+1));
        f.write(reinterpret_cast<char*>(positions),sizeof(uint32_t)*((n/samplerate)+2));
        T_bwt->save(f);
        Tr_bwt->save(f); /*Added for reverse bwt tree*/
        sampled->save(f);
        f.close();
    } else return 1;
    
    return 0;
}

FM*
FM::load(char* filename) {
    FM* newIdx = new FM();
    std::ifstream f;
    f.open(filename,std::ios::in | std::ios::binary); 
    
    if(f.is_open()) {
        info("loading FM Index from file '%s'",filename);
        f.read(reinterpret_cast<char*>(&newIdx->samplerate),sizeof(uint32_t));
        f.read(reinterpret_cast<char*>(&newIdx->sigma),sizeof(uint32_t));
        f.read(reinterpret_cast<char*>(&newIdx->I),sizeof(uint32_t));
        f.read(reinterpret_cast<char*>(&newIdx->n),sizeof(uint32_t));
        f.read(reinterpret_cast<char*>(newIdx->C),sizeof(uint32_t)*(size_uchar+1));
        f.read(reinterpret_cast<char*>(newIdx->remap),sizeof(uint8_t)*size_uchar);
        newIdx->remap_reverse = (uint8_t*) safe_malloc(sizeof(uint8_t)*(newIdx->sigma));
        f.read(reinterpret_cast<char*>(newIdx->remap_reverse),sizeof(uint8_t)*newIdx->sigma);
        newIdx->suffixes = (uint32_t*) safe_malloc(sizeof(uint32_t)*((newIdx->n/newIdx->samplerate)+1));
        f.read(reinterpret_cast<char*>(newIdx->suffixes),sizeof(uint32_t)*((newIdx->n/newIdx->samplerate)+1));
        newIdx->positions = (uint32_t*) safe_malloc(sizeof(uint32_t)*((newIdx->n/newIdx->samplerate)+2));
        f.read(reinterpret_cast<char*>(newIdx->positions),sizeof(uint32_t)*((newIdx->n/newIdx->samplerate)+2));
        newIdx->T_bwt = WaveletTreeNoptrs::load(f);
        newIdx->Tr_bwt = WaveletTreeNoptrs::load(f); /*Added the reverse bwt tree loading*/
        newIdx->sampled = BitSequenceRRR::load(f);
        f.close();
        info("samplerate = %d",newIdx->samplerate);
        info("sigma = %d",newIdx->sigma);
        info("I = %d",newIdx->I);
        info("n = %d",newIdx->n);
    } else {
        delete newIdx;
        return NULL;
    }

    return newIdx;
}

/*Count how many instances of pattern of length m appear in T*/
uint32_t
FM::count(uint8_t* pattern,uint32_t m) {
    uint8_t c = remap[pattern[m-1]]; /* map pattern to our alphabet */
    uint32_t i=m-1;
    uint32_t j = 1;
    /*sp and ep is the suffix range of the pattern*/
    uint32_t sp = C[c]; /* starting range in M from p[m-1] */
    uint32_t ep = C[c+1]-1;
	/* while there are possible occs and pattern not done */
    while (sp<=ep && i>=1) { 
      c = remap[pattern[--i]]; /* map pattern to our alphabet */
      sp = C[c] + T_bwt->rank(c, sp-1); /* LF Mapping */
      ep = C[c] + T_bwt->rank(c, ep)-1; /* LF Mapping */
      j++;
  }

  if (sp<=ep) {
      return ep-sp+1;
    } else { /*There is no appearance of P in T*/
      return 0;
  }
}

/*Given a substring P, returns an array of indexes in which P appears from that index in the string T including overlapping instances*/
uint32_t*
FM::locate(uint8_t* pattern,uint32_t m,uint32_t* matches) {
    uint32_t* locations;
    uint8_t c =  remap[pattern[m-1]];
    uint32_t i=m-1;
    
    /* count occs */
    uint32_t sp = C[c];
    uint32_t ep = C[c+1]-1;
    while (sp<=ep && i>=1) {
      c =  remap[pattern[--i]];
      sp = C[c] + T_bwt->rank(c, sp-1);
      ep = C[c] + T_bwt->rank(c, ep)-1;
  }

  if (sp<=ep) {
        /* determine positions */
        *matches = ep-sp+1; /*Number of instances of pattern in T, calculated by the range between sp and ep*/
        uint32_t locate=0;
        locations= (uint32_t*) safe_malloc((*matches)*sizeof(uint32_t));
        i=sp;
        int32_t j,dist,rank;
        while (i<=ep) {
            j=i,dist=0;
            while (!sampled->access(j)) { /*"sampled" is the Bit sequence of the wavelet tree, (maybe will be 1 if the letter in index j appears in the suffix we are looking for)*/
                c = T_bwt->access(j);
                rank = T_bwt->rank(c,j)-1;
                j = C[c]+rank; // LF-mapping
                ++dist;
            }
            locations[locate]=suffixes[sampled->rank1(j)-1]+dist;
            locate++;
            ++i;
        }
        /* locations are in SA order */
        std::sort(locations,locations+(*matches));
        return locations;
    } else {
        /* no matches */
        *matches = 0;
        return NULL;
    }

    return locations;
}

uint8_t*
FM::extract(uint32_t start,uint32_t stop)
{
    uint8_t* T;
    uint32_t m,j,skip,todo,dist;
    uint8_t c;

	/* last text pos is n-2 */
    if(stop > (this->n-1) ) stop = n-2; 
    if(start > stop) {
      return NULL;
    }

    m = stop-start+1; /* snippet len */
    T = (uint8_t*) safe_malloc( (m+1) * sizeof(uint8_t)  );

	/* determine start pos of backwards search */
     j = positions[(stop/samplerate)+1];

	/* determine distance from start pos to the text snippet we want */
    if ((stop/samplerate+1) == ((n-1)/samplerate+1)) 
        skip = n-2 - stop;
    else 
        skip = (samplerate-stop)%samplerate-1;

	/* start the backwards search */
    todo = m;
    dist = 0;
    while(todo>0) {
        c = T_bwt->access(j);
        j = C[c] + T_bwt->rank(c,j)-1;

		/* check if we are at the snippet */
        if(dist>=skip) {
            c = remap_reverse[c];
            T[todo-1] = c;
            todo--;
        }
    dist++;
    }

	/* terminate */
    T[m] = 0;

    return T;
}

uint8_t*
FM::reconstructText(uint32_t* size)
{
    uint8_t* T;
    uint8_t c;
    uint32_t j,i;
    
    T = (uint8_t*) safe_malloc( n * sizeof(uint8_t)  );
    
    j = I; /* I is sa[I] = 0 -> last sym in T */
    for(i=0;i<n;i++) {
        c = T_bwt->access(j); /* L[j] = c */
        T[n-i-1] = remap_reverse[c]; /* undo sym mapping */
        j = C[c]+T_bwt->rank(c,j)-1; /* LF-mapping: j = LF[j] */
    }
    
    if(T[n-1] == 0) *size = n-1;
    else *size = n;
    
    return T;
}

/***************************************************************************
****************************ADDED FUNCTIONS*********************************
***************************************************************************/
/*returns a vector of suffix ranges (SR) of all the appearances of the pattern in T with 1 error*/
/*TODO*/
std::vector<SR2>
FM::locate1err(uint8_t* pattern,uint32_t m) {
    uint32_t x = ceil(m / 2); /*Middle of pattern index*/
    std::vector<SR2> ans; /*Vector of ranges containing the 1 error pattern*/
    SR sr, srR; /*suffix range for SR and SR reverse*/
    SR2 sr2; /*Contains both SR and SR reverse*/
    uint32_t e; /*The character representing the error character*/
    uint32_t c;

    SR2 tmp;
    SR2 tmp_ptr;
    c = remap[pattern[m-1]]; /*The last remapped character of P*/
    sr.setRange(C[c], C[c+1] - 1); /*The initial range for the last character*/
    srR.setRange(C[c], C[c+1] - 1); /*The initial range for the last character, identical to the SR*/
    sr2.setSR(sr);
    sr2.setSRr(srR);
    /*Case 1: The error is in the first half -> P[x+1...m-1] appears in T and somewhere in P[1...x] there is 1 error*/
    tmp_ptr = StartBackwardSearch(pattern, sr2, x + 1, m - 2); /*Return the suffix range of the 2nd half of pattern (P[x+1...m-1])*/
    if (tmp_ptr.valid != 0) { /*If the 2nd part of P appears in T*/
        SR2 baseRange = tmp_ptr; /*Will hold tmp_ptr value in case BackwardSearch will return NULL so the value won't be lost*/
        for (int32_t i = x; i >= 0 ; --i) { /*If i = 0 and we reduce 1 from i, the condition i <= m-1 will no longer apply as i is an unsigned int*/
            tmp_ptr = StartBackwardSearch(pattern, baseRange, i + 1, x); /*Finding the range for P[i+1...x] assuming P[x+1...m-1] exists in T, and i is the error char*/
            if (tmp_ptr.valid == 0) continue; /*The P[i+1..x] doesn't appear in T -> Go to next iteration*/
            SR2 tmp_ptr_copy = tmp_ptr; /*Another SR2 object is used to prevent the deletion in getNewSR2 of the current found range P[0...i-1]*/
            for (e = 1; e < sigma ; e++) { /*Go over all the characters in Sigma that are not c for cP[i+1...x]*/
                if (e == remap[pattern[i]]) continue; /*e = c for cP -> go to next iteration*/
                tmp = getNewSR2(tmp_ptr_copy, e, T_bwt); /*Find the suffix range for eP[i+1...m-1]*/
                if (!tmp.isLegalRange()) continue; 
                SR2 sr_1error = tmp;
                if (i == 0)
                    ans.push_back(sr_1error); /*The error is in index 0 of P so |eP[1...m-1]| = |P| -> this pattern fits with 1 error*/
                else {
                    tmp_ptr = StartBackwardSearch(pattern, sr_1error, 0, i - 1); /*eP[i+1...m-1] exists, so we'll search for P[0..i-1] to get P[0..i-1]eP[i+1..m-1]*/
                    if (tmp_ptr.valid == 0) continue; 
                    ans.push_back(tmp_ptr); /*P[0..i-1]eP[i+1..m-1] appears in T -> we will add to the vector of 1 error matching*/
                }
            }
        }
    }

    /*Case 2: The error is in the 2nd half -> P[0..x] appears in T and somewhere in P[x+1..m-1] there is 1 error*/
    c = remap[pattern[0]]; /*The 1st remapped character of P*/
    sr.setRange(C[c], C[c+1] - 1); /*The initial range for the last character*/
    srR.setRange(C[c], C[c+1] - 1); /*The initial range for the last character, identical to the SR*/
    sr2.setSR(sr);
    sr2.setSRr(srR);
    tmp_ptr = StartForwardSearch(pattern, sr2, 1, x);/*Search a match for P[1...x] in T (we checked the ranged for P[0] above*/
    /*Note: when the line has a similar action as case one, we discarded the comment*/
    if (tmp_ptr.valid != 0) {
        SR2 baseRange = tmp_ptr;
        for (uint32_t i = x + 1 ; i <= m - 1 ; ++i) {
            tmp_ptr = StartForwardSearch(pattern, baseRange, x + 1, i - 1); /*Assuming P[0...x] exists in T -> search for P[x+1...i-1] in T*/
            if (tmp_ptr.valid == 0) continue; /*P[x+1...i-1] is not in T*/
            SR2 tmp_ptr_copy = tmp_ptr;
            for (e = 1; e < sigma; ++e) {
                if (e == remap[pattern[i]]) continue;
                tmp = getNewSR2(reverseSR_SRr(tmp_ptr_copy), e, Tr_bwt); /*Sending the reverse the order of SR and SR reverse as the forward seach is using getNewSR2 in reverse order*/
                if (!tmp.isLegalRange()) continue;
                SR2 sr_1error = reverseSR_SRr(tmp); /*Return the order of SR and SR reverse to the original state (as we switched their order for the forward search)*/
                if (i == m - 1)
                    ans.push_back(sr_1error);
                else {
                    tmp_ptr = StartForwardSearch(pattern, sr_1error, i + 1, m - 1); /*After finding a match for P[0...x..i-1]e in T, we will search for P[i+1...m-1] in T*/
                    if (tmp_ptr.valid == 0) continue;
                    ans.push_back(tmp_ptr); /*We found a match for P[0...i-1]eP[i+1..m-1] -> will be added to the vector of 1 error matching*/
                }
            }
        }
    }
    return ans;
}

/*returns a vector of suffix ranges (SR) of all the appearances of the pattern in T with 2 errors*/
std::vector<SR2>
FM::locate2err(uint8_t* pattern,uint32_t m) {
    uint32_t s1 = floor(m / 3.0); /*Middle of pattern index*/
    uint32_t s2 = m - s1;
    int32_t i, j;
    std::vector<SR2> ans; /*Vector of ranges containing the 1 error pattern*/
    SR sr, srR;
    SR2 sr2;
    SR2 tmp_ptr;
    uint32_t e2,e1; /*The character representing the error character*/

    //case a: The mismatches occur in the first 2 parts.
    uint32_t c = remap[pattern[m-1]]; /*The last remapped character of P*/
    sr.setRange(C[c], C[c+1] - 1); /*The initial range for the last character*/
    srR.setRange(C[c], C[c+1] - 1); /*The initial range for the last character, identical to the SR*/
    sr2.setSR(sr);
    sr2.setSRr(srR);

    tmp_ptr = StartBackwardSearch(pattern, sr2, s2, m - 2); /*Return the suffix range of the 3rd part of the pattern (P[x+1...m-1])*/
    if (tmp_ptr.valid != 0) {
        SR2 last_part_suffix_range = tmp_ptr;
        for (i = s2 - 1; i >= 1 /*&& i <= m - 1*/ ; --i){ //We check each sub-pattern P'[i+1...s2-1] that is correct until we find the 1st error
            tmp_ptr = StartBackwardSearch(pattern, last_part_suffix_range, i + 1, s2 - 1);
            if (tmp_ptr.valid == 0) continue; //There is no matching interval without error for the current indexes.
            SR2 last_correct_suffix_range = tmp_ptr;

            for (e2 = 1; e2 < sigma ; e2++) {
                if (e2 == remap[pattern[i]]) continue; //The current sub-pattern e2P[i+1...s2-1] fits original pattern but we are looking for an error character.
                sr2 = getNewSR2(last_correct_suffix_range, e2, T_bwt); /*Find the suffix ranges for e2P[i+1...s2-1]*/
                if (sr2.isLegalRange() == 0) continue; //Suffix range is illegal, so is irrelevant for the search.
                SR2 suffix_range_1_error = sr2;
                for (j = i - 1 ; j >= 0 ; j--) {
                    tmp_ptr = StartBackwardSearch(pattern, suffix_range_1_error, j + 1, i - 1); //continue backward search for sub-pattern P[j+1...i-1] without errors (completing the eP[i+1...m-1])
                    if (tmp_ptr.valid == 0) continue;
                    SR2 suffix_range_1_error_down_to_j = tmp_ptr;
                    for (e1 = 1 ; e1 < sigma ; e1++) {
                        if (e1 == remap[pattern[j]]) continue;
                        sr2 = getNewSR2(suffix_range_1_error_down_to_j, e1, T_bwt);
                        if (sr2.isLegalRange() == 0) continue;
                        SR2 suffix_range_2_errors = sr2;
                        if (j == 0){ //the error e1 is the 1st character, so we don't need to continue the search
                            ans.push_back(suffix_range_2_errors);
                            continue;
                        }
                        tmp_ptr = StartBackwardSearch(pattern, suffix_range_2_errors, 0, j - 1);
                        if (tmp_ptr.valid == 0) continue;
                        ans.push_back(tmp_ptr);
                    }
                }
            }
        }
    }

    //case b: Both mismatches appear in the 3rd part.
    c = remap[pattern[0]]; /*The 1st remapped character of P*/
    sr.setRange(C[c], C[c+1] - 1); /*The initial range for the 1st character*/
    srR.setRange(C[c], C[c+1] - 1); /*The initial range for the 1st character, identical to the SR*/
    sr2.setSR(sr);
    sr2.setSRr(srR);
    tmp_ptr = StartForwardSearch(pattern, sr2, 1, s2 - 1);
    if (tmp_ptr.valid != 0){ //check that the sub-pattern P[1...s2-1] is correct.
        SR2 first_two_parts_suffix_range = tmp_ptr;
        for (i = s2 ; i <= m - 2 ; i++) {
            tmp_ptr = StartForwardSearch(pattern, first_two_parts_suffix_range, s2, i - 1);
            if (tmp_ptr.valid == 0) continue;
            SR2 last_correct_suffix_range = tmp_ptr;
            for (e1 = 1 ; e1 < sigma ; e1++) {
                if (e1 == remap[pattern[i]]) continue;
                sr2 = getNewSR2(reverseSR_SRr(last_correct_suffix_range), e1, Tr_bwt);
                if (sr2.isLegalRange() == 0) continue;
                SR2 suffix_range_1_error = reverseSR_SRr(sr2);

                //Search for the 2nd error:
                for (j = i + 1 ; j <= m - 1 ; j++) {
                    tmp_ptr = StartForwardSearch(pattern, suffix_range_1_error, i + 1, j - 1);
                    if (tmp_ptr.valid == 0) continue;
                    SR2 suffix_range_1_error_up_to_j = tmp_ptr;

                    //search for the possible e2
                    for (e2 = 1 ; e2 < sigma ; e2++){
                        if (e2 == remap[pattern[j]]) continue;
                        sr2 = getNewSR2(reverseSR_SRr(suffix_range_1_error_up_to_j), e2, Tr_bwt);
                        if (sr2.isLegalRange() == 0) continue;
                        SR2 suffix_range_2_errors = reverseSR_SRr(sr2);
                        tmp_ptr = StartForwardSearch(pattern, suffix_range_2_errors, j + 1, m - 1);
                        if (tmp_ptr.valid == 0) continue;
                        ans.push_back(tmp_ptr);
                    }
                }
            }
        }
    }

    //case c: The mismatches occur in the 2nd and 3rd part.
    c = remap[pattern[0]]; /*The 1st remapped character of P*/
    sr.setRange(C[c], C[c+1] - 1); /*The initial range for the 1st character*/
    srR.setRange(C[c], C[c+1] - 1); /*The initial range for the 1st character, identical to the SR*/
    sr2.setSR(sr);
    sr2.setSRr(srR);

    tmp_ptr = StartForwardSearch(pattern, sr2, 1, s1 - 1);
    if (tmp_ptr.valid != 0) {
        SR2 first_part_suffix_range = tmp_ptr;
        for (i = s1 ; i <= s2 - 1 ; i++) {
            tmp_ptr = StartForwardSearch(pattern, first_part_suffix_range, s1, i - 1);
            if (tmp_ptr.valid == 0) continue;
            SR2 last_correct_suffix_range = tmp_ptr;
            for (e1 = 1 ; e1 < sigma ; e1++) {
                if (e1 == remap[pattern[i]]) continue;
                sr2 = getNewSR2(reverseSR_SRr(last_correct_suffix_range), e1, Tr_bwt);
                if (sr2.isLegalRange() == 0) continue;
                SR2 suffix_range_1_error = reverseSR_SRr(sr2);
                tmp_ptr = StartForwardSearch(pattern, suffix_range_1_error, i + 1, s2 - 1); //complete the search for the 2nd part.
                if (tmp_ptr.valid == 0) continue;
                SR2 suffix_range_1_error_left_parts = tmp_ptr;

                //Search for the 2nd error in the 3rd part
                for (j = s2 ; j < m ; j++) {
                    tmp_ptr = StartForwardSearch(pattern, suffix_range_1_error_left_parts, s2, j - 1);
                    if (tmp_ptr.valid == 0) continue;
                    SR2 suffix_range_1_error_up_to_j = tmp_ptr;

                    //search for the possible e2
                    for (e2 = 1 ; e2 < sigma ; e2++){
                        if (e2 == remap[pattern[j]]) continue;
                        sr2 = getNewSR2(reverseSR_SRr(suffix_range_1_error_up_to_j), e2, Tr_bwt);
                        if (sr2.isLegalRange() == 0) continue;
                        SR2 suffix_range_2_errors = reverseSR_SRr(sr2);
                        tmp_ptr = StartForwardSearch(pattern, suffix_range_2_errors, j + 1, m - 1);
                        if (tmp_ptr.valid == 0) continue;
                        ans.push_back(tmp_ptr);
                    }
                }
            }
        }
    }


    //case d: The mismtaches occur in the 1st and 3rd part.
    c = remap[pattern[s1]]; /*The remapped character of the 2nd part of P*/
    sr.setRange(C[c], C[c+1] - 1); /*The initial range for the character*/
    srR.setRange(C[c], C[c+1] - 1); /*The initial range for the character, identical to the SR*/
    sr2.setSR(sr);
    sr2.setSRr(srR);

    tmp_ptr = StartForwardSearch(pattern, sr2, s1 + 1, s2 - 1);
    if (tmp_ptr.valid != 0) {
        SR2 middle_part_range = tmp_ptr;
        for(i = s1 - 1 ; i >= 0 ; i--){
            tmp_ptr = StartBackwardSearch(pattern, middle_part_range, i + 1, s1 - 1);
            if (tmp_ptr.valid == 0) continue;
            SR2 suffix_back_no_error = tmp_ptr;

            for(e1 = 1 ; e1 < sigma ; e1++){
                if(e1 == remap[pattern[i]]) continue;
                sr2 = getNewSR2(suffix_back_no_error,e1,T_bwt);
                if (sr2.isLegalRange() == 0) continue;
                SR2 suffix_back_1_error = sr2;
                SR2 full_left_parts_1_error;
                if (i != 0) {
                    tmp_ptr = StartBackwardSearch(pattern, suffix_back_1_error, 0, i - 1);
                    if (tmp_ptr.valid == 0) continue;
                    full_left_parts_1_error = tmp_ptr;
                }
                else
                    full_left_parts_1_error = suffix_back_1_error;

                for (j = s2 ; j < m ; j++) {
                    tmp_ptr = StartForwardSearch(pattern, full_left_parts_1_error, s2, j- 1);
                    if (tmp_ptr.valid == 0) continue;
                    SR2 third_part_no_error = tmp_ptr;

                    for (e2 = 1 ; e2 < sigma ; e2++) {
                        if (e2 == remap[pattern[j]]) continue;
                        sr2 = getNewSR2(reverseSR_SRr(third_part_no_error), e2, Tr_bwt);
                        if (sr2.isLegalRange() == 0) continue;
                        SR2 third_part_1_error = reverseSR_SRr(sr2);
                        tmp_ptr = StartForwardSearch(pattern, third_part_1_error, j + 1, m - 1);
                        if (tmp_ptr.valid == 0) continue;
                        ans.push_back(tmp_ptr);
                    }
                }
            }
        }
    }
    
    return ans;
}

/*Given a suffix range of matching pattern, the function will return the starting index of each match and change matches variable value to the number of matches*/
uint32_t*
FM::getIndexes(SR sr, uint32_t *matches) {
    *matches = sr.getEp() - sr.getSp() + 1;
    uint32_t *ans = (uint32_t*) safe_malloc((*matches)*sizeof(uint32_t));
    uint32_t counter = 0, diff, curr;
    for (uint32_t i = sr.getSp() ; i <= sr.getEp() ; i++){
        diff = 0;
        curr = i;
        while (sampled->access(curr) == 0) {
            curr = C[T_bwt->access(curr)] + T_bwt->rank(T_bwt->access(curr),curr-1); /*slides 22-24, finding LF if the current index is not sampled*/
            diff++;
        }
        ans[counter++] = suffixes[sampled->rank1(curr)-1] + diff; /*The current index is the index of the previews sampled suffix + the difference we found*/
    }
    return ans;
}


SR2
FM::reverseSR_SRr(SR2 sr2){
    sr2 = SR2(sr2.getSRr(), sr2.getSR());
    sr2.valid = 1;
    return sr2;
}

/*Return the number of suffixes Pd, s.t. d < c*/
size_t FM::smallerCharSuffixes(uint8_t c, uint32_t bound, WaveletTreeNoptrs *T){
    if (bound == 0 || bound == UINT32_MAX)
        return 0;
    size_t ans = 0;
    for (uint8_t i = 0 ; i < c ; i++) {
        ans += T->rank(i, bound); /*Uses the rank function of the Wavelet tree for each i < c*/
    }
    return ans;
}


SR2
FM::StartBackwardSearch(uint8_t* pattern, SR2 sr2, int32_t sp, int32_t ep) {
    if (sp > ep){ /*If the range is illegal we will return the SR and SRr structure as is*/
        sr2.valid = 1;
        return sr2;
    }
    for (int32_t i = ep ; 
        sr2.getSR().getSp() <= sr2.getSR().getEp() && /*The range for the given SR is legal*/
        sr2.getSRr().getSp() <= sr2.getSRr().getEp() && /*The range for the given SRr is legal*/
        i >= sp &&
        i <= ep; /*i is in the range between sp and ep*/
        --i) {
            sr2 = getNewSR2(sr2, remap[pattern[i]], T_bwt); /*Get the range for cP where the index of c is i*/
    }
    if (!sr2.isLegalRange()){ /*If the sp and ep switched then the pattern does not appear in the string*/
        sr2.valid = 0;
        return sr2;
    }
    sr2.valid = 1;
    return sr2; /*return the new interval of SR and SRr containing the final range of the pattern*/
}


/*Return the next SR and SRr given a current range and the reverse range*/
/*Assume [s',e'] is the next range of [s,e] for SR
  Assume [sR',eR'] is the next range of [sR,eR] for SRr*/
size_t
FM::getRank(WaveletTreeNoptrs *T, uint8_t c, int32_t range){
    if (range <= 0) return 0;
    return T->rank(c,range);
}

SR2
FM::getNewSR2(SR2 sr2, uint8_t c, WaveletTreeNoptrs *T) {
    /*Simplifying the use of SR and SRr*/
    SR srR_interval = sr2.getSRr();
    SR sr_interval = sr2.getSR();
    /*sR' = sR + SIGMA(d < c):{{ rank(d,Tr_bwt,e) - rank(d,Tr_bwt,s-1), slide 33*/
    srR_interval.setSp(srR_interval.getSp() + smallerCharSuffixes(c,sr_interval.getEp(),T) -  smallerCharSuffixes(c,sr_interval.getSp() - 1, T));
    /*eR' = sR' + (e' - s'), slide 33*/
    srR_interval.setEp(/*new SP*/srR_interval.getSp() + T->rank(c,sr_interval.getEp()) - /*getRank(T,c,sr_interval.getSp()-1)*/ T->rank(c,sr_interval.getSp() - 1) - 1);
    /*Get the next [s',e'] given the range [s,e] using the regular backward seach method for the next substring cP*/
    sr_interval = getNewBackSR(sr_interval, c, T);
    SR2 ans = SR2(sr_interval.getSp(), sr_interval.getEp(), srR_interval.getSp(), srR_interval.getEp());
    ans.valid = 1;
    return ans;
}


SR
FM::getNewBackSR(SR sr, uint8_t c, WaveletTreeNoptrs *T) {
    /*
    Given range [l,r] and character c
    l' = C[c] + rank(c,T,l-1)
    r' = C[c] + rank(c,T,r)
    slide 14 (note: the slide shows + 1 on l', we don't use it at the moment)
    */
    sr.setSp(C[c] +  T->rank(c,sr.getSp()-1));
    sr.setEp(C[c] +  T->rank(c,sr.getEp())-1);
    return sr;
}


SR2
FM::StartForwardSearch(uint8_t* pattern, SR2 sr2, int32_t sp, int32_t ep) { /*given the su*/
    if (sp > ep){ /*If the range is illegal we will return the SR and SRr structure as is*/
        sr2.valid = 1;
        return sr2;
    }
    sr2 = reverseSR_SRr(sr2);
    for (uint32_t i = sp ; 
        sr2.getSR().getSp() <= sr2.getSR().getEp() && /*The range for the given SR is legal*/
        sr2.getSRr().getSp() <= sr2.getSRr().getEp() && /*The range for the given SRr is legal*/
        i >= sp &&
        i <= ep; /*i is in the range between sp and ep*/
        ++i) {
            sr2 = getNewSR2(sr2, remap[pattern[i]], Tr_bwt); /*Get the range for Pc where the index of c is i*/
        }
    if (!sr2.isLegalRange()){ /*If the sp and ep switched then the pattern does not appear in the string*/
        sr2.valid = 0;
        return sr2;
    }
    sr2.valid = 1;
    return reverseSR_SRr(sr2); /*return the new interval of SR and SRr containing the final range of the pattern*/
}


/*Creates the reverse remapped string Xr for a given X string*/
uint8_t*
FM::reverse(uint8_t* X, uint32_t n) {
    uint32_t i;
    uint8_t* Xr = (uint8_t*) malloc (sizeof(uint8_t) * n);
    for (i = 0 ; i < n - 1 ; i++) /*Not including the \0 at the end of X*/
        Xr[i] = X[n-i-2]; /*minus 2 for the '\0' and the index beginning at 0 and ends at n-1*/
    Xr[i] = '\0';
    return Xr;
}
