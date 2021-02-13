
# include <fstream>
# include <stdlib.h>
# include <string>
# include <sstream>
# include <iostream>
# include <unordered_map>
# include <unordered_set>
# include <vector>
# include <math.h>
# include <algorithm>
# include <iomanip>
# include <boost/math/distributions/chi_squared.hpp>



// Global variables for various functions.

int buffer = 10;
int globorder = 2;
int knum = 4;
int minseglen = 64;


// Cluster struct to maintain distribution counts, total length, and number of segments (normalization factor)

struct cluster {
    float length;
    float segcount;
    std::vector<std::vector<float>> kmers;
    std::string nucleo;

};

// Distribution for dMAX equation.

boost::math::chi_squared chi2dist(48);
boost::math::chi_squared alt_chi2dist(63);

// Distribution library for segments.

std::unordered_map<std::string,std::vector<std::vector<float>>> distroz;

// Simple Indices for converting kmers.

std::unordered_map<std::string,int> indexmer {
    {"AA",0},
    {"AT",1},
    {"AG",2},
    {"AC",3},
    {"TA",4},
    {"TT",5},
    {"TG",6},
    {"TC",7},
    {"GA",8},
    {"GT",9},
    {"GG",10},
    {"GC",11},
    {"CA",12},
    {"CT",13},
    {"CG",14},
    {"CC",15},
};

std::unordered_map<int,std::string> rindexmer {
    {0,"AA"},
    {1,"AT"},
    {2,"AG"},
    {3,"AC"},
    {4,"TA"},
    {5,"TT"},
    {6,"TG"},
    {7,"TC"},
    {8,"GA"},
    {9,"GT"},
    {10,"GG"},
    {11,"GC"},
    {12,"CA"},
    {13,"CT"},
    {14,"CG"},
    {15,"CC"},
};

// Nucleotide to position translator for distributions.

std::unordered_map<char,int> nucloc = {
    {'A',1},
    {'T',2},
    {'G',3},
    {'C',4}
};

std::unordered_map<int,char> rnucloc = {
    {1,'A'},
    {2,'T'},
    {3,'G'},
    {4,'C'}
};

// Zero-safe log function.


double safelog(float count1,float count2) {
    double f;
    double x,y;
    y = count2;
    x = count1;
    if (y == 0.0) {
        f = 0.0;
    }
    else if ((x/y) == 1.0) {
        f = 0.0;
    }
    else if ((x/y) == 0.0) {
        f = 0.0;
    }
    else {
        f = ((x/y) * log2(x/y));
    }
    
    return f;
}

// IslandCAFE Parameters for dmax

double dmax(double mjsd,float glen) {
    double dmjsd;
    double neff = neff = -7.66 + 2.39*(log(glen));
    double beta = beta = 0.841 + 0.0029*(log(glen));
    double interior = (2*glen)*log(2)*mjsd*beta;
    if (mjsd > 0) {
        dmjsd = pow(boost::math::cdf(chi2dist,interior),neff);
    }
    else {
        dmjsd = 1.0;
    }
    return dmjsd;

}

// Alternative parameter dmax
// Easier to just rename this to dmax (and the other to alt_dmax) if you want to test their accuracy.

double alt_dmax(float mjsd,float glen) {
    double dmjsd;
    double neff = neff = -2.44732 + 1.13049*(log(glen));
    double beta = beta = 1.025173 + 0.002304*(log(glen));
    double interior = (2*glen)*log(2)*mjsd*beta;
    if (mjsd > 0) {
        dmjsd = pow(boost::math::cdf(alt_chi2dist,interior),neff);
    }
    else {
        dmjsd = 1.0;
    }
    return dmjsd;

}

// Fasta import and filter function.
// Skips id lines, returns only standard ATGC nucleotides in uppercase.

std::string fasta(std::string infile) {
    std::ifstream fna(infile);
    std::string line;
    std::string genome;
    while(std::getline(fna,line)) {
        if (line[0] == '>') {
            continue;
        }
        for(char const &c: line) {
            if (c == 'A' || c == 'a') {genome+=('A');}
            if (c == 'T' || c == 't') {genome+=('T');}
            if (c == 'G' || c == 'g') {genome+=('G');}
            if (c == 'C' || c == 'c') {genome+=('C');}
        }
        
        
    }
    return genome;
}

// Iterative mJSD function through reverse chain algorithm.

std::pair<double,int> mJSDi(std::string gstring) {
    double BMJSD,MJSD;
    float slength = gstring.length();
    float lenA = buffer;
    float lenB = slength-lenA;
    float bestloc = lenA;
    double HA = 0;
    double HB = 0;
    double HG = 0;


    
    // Init count structs.
    // Store as kmer <total,a,t,g,c>
    
    std::vector<std::vector<float>> distA;
    std::vector<std::vector<float>> distB;
    std::vector<std::vector<float>> distG;
    std::vector<std::vector<float>> bestdistA;
    std::vector<std::vector<float>> bestdistB;
    distA.resize(16,std::vector<float>(5));
    distB.resize(16,std::vector<float>(5));
    distG.resize(16,std::vector<float>(5));
    
    // Add initial counts from first substring to A distro.
    // Add same counts to combined substring G distro.

    for(int i=0;i<(buffer-globorder);i++) {
        std::string&& kmer = gstring.substr(i,2);
        char afterk = gstring[i+globorder];
        int pos;
        if (afterk == 'A') pos=1;
        if (afterk == 'T') pos=2;
        if (afterk == 'G') pos=3;
        if (afterk == 'C') pos=4;
        distA[indexmer[kmer]][pos]++;
        distA[indexmer[kmer]][0]++;
        distG[indexmer[kmer]][pos]++;
        distG[indexmer[kmer]][0]++;
    }
    
    // Add initial counts from first substring to B distro.
    // Add same counts to combined substring G distro.    

    for(int i=buffer;i<(slength-globorder);i++) {
        std::string&& kmer = gstring.substr(i,2);
        char afterk = gstring[i+globorder];
        int pos;
        if (afterk == 'A') pos=1;
        if (afterk == 'T') pos=2;
        if (afterk == 'G') pos=3;
        if (afterk == 'C') pos=4;
        distB[indexmer[kmer]][pos]++;
        distB[indexmer[kmer]][0]++;
        distG[indexmer[kmer]][pos]++;
        distG[indexmer[kmer]][0]++;
    }    

    // Calculate initial entropy.
    
    for(int ix=0;ix<=15;ix++) {
        double HA1 = 0;
        double HB1 = 0;
        double HG1 = 0;
        for(int i=1;i<=knum;i++) {
            HG1 += safelog(distG[ix][i],distG[ix][0]);
            HA1 += safelog(distA[ix][i],distA[ix][0]);
            HB1 += safelog(distB[ix][i],distB[ix][0]);
        }
        HG += ((distG[ix][0]/(slength-2)) * HG1);
        HA += ((distA[ix][0]/(lenA-2)) * HA1);
        HB += ((distB[ix][0]/(lenB-2)) * HB1);
        
    }
    
    float wa = lenA/slength;
    float wb = lenB/slength;
    HA = -HA;
    HB = -HB;
    HG = -HG;
    BMJSD = HG - (wa*HA) - (wb*HB);

    
    // Calculate iterative entropies.
    
    double ZA = HA;
    double ZB = HB;
    double ZG = HG;
    float iterations = lenB - buffer;
    
    for(int i=0;i<iterations;i++) {
        ZA = -ZA;
        ZB = -ZB;
        ZG = -ZG;
        double ZA1 = 0;
        double ZB1 = 0;
        double ZG1 = 0;
        double ZG2 = 0;


        std::string&& gainmer = gstring.substr(i+buffer-globorder,2);
        std::string&& lossmer = gstring.substr((i+buffer),2);
        char gnuc = gstring[i+buffer];
        char lnuc = gstring[i+buffer+globorder];
        
// When gain/loss between each substring is the same kmer, only perform the ZG adjustment once.
        
        if (gainmer == lossmer) {
            for(int i2=1;i2<=knum;i2++) {
                ZA1 += safelog(distA[indexmer[gainmer]][i2],distA[indexmer[gainmer]][0]);
                ZB1 += safelog(distB[indexmer[lossmer]][i2],distB[indexmer[lossmer]][0]);
                ZG1 += safelog(distG[indexmer[lossmer]][i2],distG[indexmer[lossmer]][0]);
                ZG2 += safelog(distG[indexmer[gainmer]][i2],distG[indexmer[gainmer]][0]);
            }
            ZA -= ((distA[indexmer[gainmer]][0]/(lenA-2)) * ZA1);
            ZB -= ((distB[indexmer[lossmer]][0]/(lenB-2)) * ZB1);
            ZG -= ((distG[indexmer[lossmer]][0]/(slength-2)) * ZG1);
        }
        else {
            for(int i2=1;i2<=knum;i2++) {
                ZA1 += safelog(distA[indexmer[gainmer]][i2],distA[indexmer[gainmer]][0]);
                ZB1 += safelog(distB[indexmer[lossmer]][i2],distB[indexmer[lossmer]][0]);
                ZG1 += safelog(distG[indexmer[lossmer]][i2],distG[indexmer[lossmer]][0]);
                ZG2 += safelog(distG[indexmer[gainmer]][i2],distG[indexmer[gainmer]][0]);
            }
            ZA -= ((distA[indexmer[gainmer]][0]/(lenA-2)) * ZA1);
            ZB -= ((distB[indexmer[lossmer]][0]/(lenB-2)) * ZB1);
            ZG -= ((distG[indexmer[lossmer]][0]/(slength-2)) * ZG1);
            ZG -= ((distG[indexmer[gainmer]][0]/(slength-2)) * ZG2);
        }

        
        ZA = (ZA * (lenA-2));
        ZB = (ZB * (lenB-2));
        lenA++;
        lenB--;
        ZA = (ZA / (lenA-2));
        ZB = (ZB / (lenB-2));
        
        
        distA[indexmer[gainmer]][0]++;
        distB[indexmer[lossmer]][0]--;
        distA[indexmer[gainmer]][nucloc[gnuc]]++;
        distB[indexmer[lossmer]][nucloc[lnuc]]--;
        
        distG[indexmer[gainmer]][0]++;
        distG[indexmer[lossmer]][0]--;
        distG[indexmer[gainmer]][nucloc[gnuc]]++;
        distG[indexmer[lossmer]][nucloc[lnuc]]--;
        
        ZA1 = 0;
        ZB1 = 0;
        ZG1 = 0;
        ZG2 = 0;



        if (gainmer == lossmer) {
            for(int i3=1;i3<=knum;i3++) {
                ZA1 += safelog(distA[indexmer[gainmer]][i3],distA[indexmer[gainmer]][0]);
                ZB1 += safelog(distB[indexmer[lossmer]][i3],distB[indexmer[lossmer]][0]);
                ZG1 += safelog(distG[indexmer[lossmer]][i3],distG[indexmer[lossmer]][0]);
                ZG2 += safelog(distG[indexmer[gainmer]][i3],distG[indexmer[gainmer]][0]);
            }
            ZA += ((distA[indexmer[gainmer]][0]/(lenA-2)) * ZA1);
            ZB += ((distB[indexmer[lossmer]][0]/(lenB-2)) * ZB1);
            ZG += ((distG[indexmer[lossmer]][0]/(slength-2)) * ZG1);
        }
        else {
            for(int i3=1;i3<=knum;i3++) {
                ZA1 += safelog(distA[indexmer[gainmer]][i3],distA[indexmer[gainmer]][0]);
                ZB1 += safelog(distB[indexmer[lossmer]][i3],distB[indexmer[lossmer]][0]);
                ZG1 += safelog(distG[indexmer[lossmer]][i3],distG[indexmer[lossmer]][0]);
                ZG2 += safelog(distG[indexmer[gainmer]][i3],distG[indexmer[gainmer]][0]);
            }
            ZA += ((distA[indexmer[gainmer]][0]/(lenA-2)) * ZA1);
            ZB += ((distB[indexmer[lossmer]][0]/(lenB-2)) * ZB1);
            ZG += ((distG[indexmer[gainmer]][0]/(slength-2)) * ZG2);
            ZG += ((distG[indexmer[lossmer]][0]/(slength-2)) * ZG1);
        }

        wa = lenA/slength;
        wb = lenB/slength;
        ZA = -ZA;
        ZB = -ZB;
        ZG = -ZG;
        HG = 0;

        double MJSD = ZG - (wa*ZA) - (wb*ZB);




    if (MJSD > BMJSD) {
        BMJSD = MJSD;
        bestloc = lenA;
        bestdistA = distA;
        bestdistB = distB;
    }
        
        
    }

    distroz[gstring.substr(0,bestloc)] = bestdistA;
    distroz[gstring.substr(bestloc)] = bestdistB;
    BMJSD = dmax(BMJSD,slength);
    
    return std::make_pair(BMJSD,bestloc);
}

// Debugging distribution printer.

void dprinter(std::vector<std::vector<float>> dist) {
    for(int ix=0;ix<=15;ix++) {
        for(int i=1;i<=knum;i++) {
            std::cout << rindexmer[ix] << nucloc[i] << ": " << dist[ix][i] << std::endl;
        }
    }
}

// Segmentation function.
// Takes in string, recursively segments until impossible.  Returns vector of segments (in original order).

std::vector<std::string> segmentation(std::string gnome,float threshold) {
    std::vector<std::string> tempsegments = {gnome};
    std::vector<std::string> segments;    
    std::pair <double,int> segtest;
    std::string sub1,sub2;

    while (tempsegments.size() != 0) {
        auto it = std::prev(tempsegments.end());
        if ((*it).length() < minseglen) {
            #Fix Here
                    if ((*it).length() < minseglen) {
            std::vector<std::vector<float>> distMini;
            distMini.resize(16,std::vector<float>(5));
            for(int i=0;i<((*it).length()-globorder);i++) {
                std::string kmer = (*it).substr(i,2);
                char afterk = (*it)[i+globorder];
                int pos;
                if (afterk == 'A') pos=1;
                if (afterk == 'T') pos=2;
                if (afterk == 'G') pos=3;
                if (afterk == 'C') pos=4;
                distMini[indexmer[kmer]][pos]++;
                distMini[indexmer[kmer]][0]++;

            }
            distroz[(*it)] = distMini;
            #Fix Above
            segments.push_back(std::move(*it));
            tempsegments.pop_back();
            continue;
        }
        segtest = mJSDi(*it);
        if (segtest.first >= threshold) {
            sub1 = (*it).substr(0,segtest.second);
            sub2 = (*it).substr(segtest.second);
            tempsegments.pop_back();
            tempsegments.push_back(sub1);
            tempsegments.push_back(sub2);
        }
        else {
            segments.push_back(std::move(*it));
            tempsegments.pop_back();
            continue;
        }
    }

    std::reverse(segments.begin(),segments.end());
    
    return segments;
}

// Full Entropy (3x Distribution) function with cluster support.
// Only compatible with iteratively segmented clusters.

std::vector<double> entropy(cluster clusseg,cluster clusseg2) {
    std::vector<std::vector<float>> distA = clusseg.kmers;
    std::vector<std::vector<float>> distB = clusseg2.kmers;
    std::vector<std::vector<float>> distG,distA2,distB2;
    distG.resize(16,std::vector<float>(5));
    distA2.resize(16,std::vector<float>(5));
    distB2.resize(16,std::vector<float>(5));
    float norm = clusseg.segcount;
    float norm2 = clusseg2.segcount;
    float slen = clusseg.length;
    float slen2 = clusseg2.length;
    double HA = 0;
    double HB = 0;
    double HG = 0;
// Get total counts from all strings in each cluster.

    for(int x=0;x<=15;x++) {
        for(int y=0;y<=knum;y++) {
            distA2[x][y] = ((distA[x][y])/norm);
            distB2[x][y] = ((distB[x][y])/norm2);
            distG[x][y] = distA2[x][y] + distB2[x][y];
        }
    }
    
    float slenT = (slen/norm) + (slen2/norm2);

    
    for(int ix=0;ix<=15;ix++) {
        double HA1 = 0;
        double HB1 = 0;
        double HG1 = 0;
        for(int ixx=1;ixx<=knum;ixx++) {
            HA1 += safelog(distA2[ix][ixx],distA2[ix][0]);
            HB1 += safelog(distB2[ix][ixx],distB2[ix][0]);
            HG1 += safelog(distG[ix][ixx],distG[ix][0]); 
        }
        HA += (((distA2[ix][0]) / (slen/norm)) * HA1);
        HB += (((distB2[ix][0]) / (slen2/norm2)) * HB1);
        HG += (((distG[ix][0]) / slenT) * HG1);
    }
    
    return {HA,slen/norm,HB,slen2/norm2,HG,slenT};
}

// Full Entropy (3x Distribution) function with cluster support.
// Only compatible with iteratively segmented clusters.
// This alteration treats clusters as entire segments, not averaging.  Useful
// for using downstream of horizontal clustering, as these are re-merged pieces.

std::vector<double> entropysegment(cluster clusseg,cluster clusseg2) {
    std::vector<std::vector<float>> distA = clusseg.kmers;
    std::vector<std::vector<float>> distB = clusseg2.kmers;
    std::vector<std::vector<float>> distG,distA2,distB2;
    distG.resize(16,std::vector<float>(5));
    distA2.resize(16,std::vector<float>(5));
    distB2.resize(16,std::vector<float>(5));
    float norm = 1.0;
    float norm2 = 1.0;
    float slen = clusseg.length;
    float slen2 = clusseg2.length;
    double HA = 0;
    double HB = 0;
    double HG = 0;
// Get total counts from all strings in each cluster.

    for(int x=0;x<=15;x++) {
        for(int y=0;y<=knum;y++) {
            distA2[x][y] = ((distA[x][y])/norm);
            distB2[x][y] = ((distB[x][y])/norm2);
            distG[x][y] = distA2[x][y] + distB2[x][y];
        }
    }
    
    float slenT = (slen/norm) + (slen2/norm2);

    
    for(int ix=0;ix<=15;ix++) {
        double HA1 = 0;
        double HB1 = 0;
        double HG1 = 0;
        for(int ixx=1;ixx<=knum;ixx++) {
            HA1 += safelog(distA2[ix][ixx],distA2[ix][0]);
            HB1 += safelog(distB2[ix][ixx],distB2[ix][0]);
            HG1 += safelog(distG[ix][ixx],distG[ix][0]); 
        }
        HA += (((distA2[ix][0]) / (slen/norm)) * HA1);
        HB += (((distB2[ix][0]) / (slen2/norm2)) * HB1);
        HG += (((distG[ix][0]) / slenT) * HG1);
    }
    
    return {HA,slen/norm,HB,slen2/norm2,HG,slenT};
}

double dmJSDcalc(std::vector<double> entropies) {
    double HA = -(entropies[0]);
    float lenA = entropies[1];
    double HB = -(entropies[2]);
    float lenB = entropies[3];
    double HG = -(entropies[4]);
    float lenG = entropies[5];
    
    double wa = lenA/lenG;
    double wb = lenB/lenG;

    double mjsd = HG - (wa*HA) - (wb*HB);
    mjsd = dmax(mjsd,lenG);
    return mjsd;
    
}

double mJSDcalc(std::vector<double> entropies) {
    double HA = -(entropies[0]);
    float lenA = entropies[1];
    double HB = -(entropies[2]);
    float lenB = entropies[3];
    double HG = -(entropies[4]);
    float lenG = entropies[5];
    
    double wa = lenA/lenG;
    double wb = lenB/lenG;

    double mjsd = HG - (wa*HA) - (wb*HB);
    return mjsd;
    
}

cluster clustermerge(cluster clustera,cluster clusterb) {
    cluster clusterab;
    std::vector<std::vector<float>> distG;
    distG.resize(16,std::vector<float>(5));
    for(int x=0;x<=15;x++) {
        for(int y=0;y<=knum;y++) {
            distG[x][y] = (clustera.kmers[x][y]) + (clusterb.kmers[x][y]);
        }
    }
    clusterab.kmers = distG;
    clusterab.length = clustera.length + clusterb.length;
    clusterab.segcount = clustera.segcount + clusterb.segcount;
    clusterab.nucleo = clustera.nucleo + clusterb.nucleo;

    return clusterab;    
}

// Initial Positional (Adjacency) Clustering Function

std::vector<cluster> clusterzontal (std::vector<std::string> segments,float threshold) {
    std::vector<cluster> grpsegments;
    double dmjsd;
    for(int i=0;i<segments.size();++i) {
        cluster clusseg;
        clusseg.kmers = distroz[segments[i]];
        clusseg.length = segments[i].length();
        clusseg.segcount = 1.0;
        clusseg.nucleo = segments[i];
        grpsegments.push_back({clusseg});
    }
    for(int i=0;i<grpsegments.size()-1;++i) {
        dmjsd = dmJSDcalc(entropy(grpsegments[i],grpsegments[i+1]));
        if (dmjsd < threshold) {
            grpsegments[i] = clustermerge(grpsegments[i],grpsegments[i+1]);
            grpsegments.erase(grpsegments.begin() + i + 1);
            --i;
        }
    }
    return grpsegments;
}

// Secondary Positional (Adjacency) Clustering Function

std::vector<cluster> clusterzontal2 (std::vector<cluster> grpsegments,float threshold) {
    double dmjsd;
    for(int i=0;i<grpsegments.size()-1;++i) {
        dmjsd = dmJSDcalc(entropy(grpsegments[i],grpsegments[i+1]));
        if (dmjsd < threshold) {
            grpsegments[i] = clustermerge(grpsegments[i],grpsegments[i+1]);
            grpsegments.erase(grpsegments.begin() + i + 1);
            --i;
        }
    }
    return grpsegments;
}

// Tertiary Non-Adjacent Group Clustering Function

std::vector<cluster> clusterbeyond (std::vector<cluster> grpsegments,float threshold) {
    double dmjsd;
    std::string dubnucleo;
    for(int i=0;i<grpsegments.size()-1;++i) {
        for(int i2=i+1;i2<grpsegments.size();++i2) {
            dmjsd = dmJSDcalc(entropy(grpsegments[i],grpsegments[i2]));
            if (dmjsd < threshold) {
                grpsegments[i] = clustermerge(grpsegments[i],grpsegments[i2]);
                grpsegments.erase(grpsegments.begin() + i2);
                i2 = i;
            }
        }
    }
    return grpsegments;
}

// Matrix manipulator function to remove columns and "reset" rows for cluster memory.

std::vector<std::vector<double>> matrixShift (std::vector<std::vector<double>> inmatrix, int row, int col) {
    for(int i=0;i<inmatrix[row].size();++i) {
        inmatrix[row][i] = 99.0;
    }
    for(int i=0;i<inmatrix.size();++i) {
        inmatrix[i][row] = 99.0;
        inmatrix[i].erase(inmatrix[i].begin() + col);
    }
    inmatrix.erase(inmatrix.begin() + col);
    return inmatrix;
}

// Tertiary Non-Adjacent Group Clustering Function with best MJSD per round.

std::vector<cluster> clusterbeyondX (std::vector<cluster> grpsegments,float threshold) {
    double dmjsd;
    bool change = true;
    std::vector<std::vector<double>> matrix(grpsegments.size(),std::vector<double>(grpsegments.size(),99.0));
    while(change == true) {
        double bestlow = 5;
        int bestdex,bestdex2;
        change = false;
        for(int i=0;i<grpsegments.size()-1;++i) {
            for(int i2=i+1;i2<grpsegments.size();++i2) {
                if (matrix[i][i2] > 90) {
                    dmjsd = dmJSDcalc(entropy(grpsegments[i],grpsegments[i2]));
                    matrix[i][i2] = dmjsd;
                }
                else {
                    dmjsd = matrix[i][i2];
                }
                if ((dmjsd < threshold) && (dmjsd < bestlow)) {
                    bestlow = dmjsd;
                    bestdex = i;
                    bestdex2 = i2;
                }
            }
    
        }
        if (bestlow < threshold) {
            grpsegments[bestdex] = clustermerge(grpsegments[bestdex],grpsegments[bestdex2]);
            grpsegments.erase(grpsegments.begin() + bestdex2);
            matrix = matrixShift(matrix,bestdex,bestdex2);
            change = true;
        }
    }
    return grpsegments;
}

// Printing FASTA function.

void fasfilePrinter(std::string outfile,std::vector<cluster> clusters) {
    std::ofstream fasfile;
    fasfile.open(outfile);
    for(int i=0;i<clusters.size();i++){
        fasfile << ">" << outfile << "." << i << "\n";
        fasfile << clusters[i].nucleo;
        fasfile << "\n";
    }
    fasfile.close();
}

// Prints mcl compatible JSD file.

void netfilePrinter(std::string outfile,std::vector<cluster> clusters) {
    double mjsd;
    std::ofstream fasfile;
    fasfile.open(outfile);
    for(int i=0;i<clusters.size();i++){
        fasfile << ">" << outfile << "." << i << "\n";
        fasfile << clusters[i].nucleo;
        fasfile << "\n";
    }
    fasfile << "###...###\n";
    for(int i=0;i<clusters.size()-1;i++) {
        for(int i2=i+1;i2<clusters.size();i2++) {
            mjsd = dmJSDcalc(entropysegment(clusters[i],clusters[i2]));
            fasfile << mjsd << "\n";
        }
    }
    fasfile.close();
}

// Prints out distance matrix for use with downstream outsourced clustering options.

void disMatrixPrint(std::string outfile,std::string outfile2,std::vector<cluster> clusters) {
    double mjsd;
    const int m = clusters.size();
    float matrix[m][m];
    std::ofstream fasfile;
    std::ofstream matxfile;
    fasfile.open(outfile);
    matxfile.open(outfile2);
    for(int i=0;i<clusters.size();i++){
        fasfile << ">" << outfile << "." << i << "\n";
        fasfile << clusters[i].nucleo;
        fasfile << "\n";
    }
    for(int i=0;i<clusters.size()-1;i++) {
        for(int i2=i+1;i2<clusters.size();i2++) {
            mjsd = mJSDcalc(entropysegment(clusters[i],clusters[i2]));
            matrix[i][i2] = mjsd;
            matrix[i2][i] = mjsd;
            matrix[i][i] = 0;
            matrix[i2][i2] = 0;
        }
    }
    for(int i=0;i<m;i++) {
        for(int i2=0;i2<m;i2++) {
            if (i2 == m-1) {
                matxfile << matrix[i][i2] << "\n";
            }
            else {
                matxfile << matrix[i][i2] << "\t";
            }
        }
    }
}

void simMatrixPrint(std::string outfile,std::string outfile2,std::vector<cluster> clusters) {
    double mjsd;
    const int m = clusters.size();
    float matrix[m][m];
    std::ofstream fasfile;
    std::ofstream matxfile;
    fasfile.open(outfile);
    matxfile.open(outfile2);
    for(int i=0;i<clusters.size();i++){
        fasfile << ">" << outfile << "." << i << "\n";
        fasfile << clusters[i].nucleo;
        fasfile << "\n";
    }
    for(int i=0;i<clusters.size()-1;i++) {
        for(int i2=i+1;i2<clusters.size();i2++) {
            mjsd = dmJSDcalc(entropysegment(clusters[i],clusters[i2]));
            matrix[i][i2] = mjsd;
            matrix[i2][i] = mjsd;
            matrix[i][i] = 1;
            matrix[i2][i2] = 1;
        }
    }
    for(int i=0;i<m;i++) {
        for(int i2=0;i2<m;i2++) {
            if (i2 == m-1) {
                matxfile << matrix[i][i2] << "\n";
            }
            else {
                matxfile << matrix[i][i2] << "\t";
            }
        }
    }
}

// Pull in range list for looping.

void getRanges(std::string rangeFile, std::vector<std::string> & range) {
    std::ifstream in(rangeFile.c_str());
    std::string ran;
    while (std::getline(in,ran)) {
        if (ran.size() > 0) {
            range.push_back(ran);
        }
    }
    in.close();
}

// Main Function

int main(int argc, char *argv[]) {
    // Intake arguments.
    std::string infile = argv[1];
    std::string outfile;
    std::vector<std::string> range;
    getRanges("Ranges.txt",range);
    float t1 = atof(argv[2]);
    float t2 = atof(argv[3]);
    double t3;
    outfile = argv[4];
    std::cout << "Importing FASTA" << std::endl;
    std::string gnome = fasta(infile);
    std::cout << "Segmenting Genome" << std::endl;
    std::vector<std::string> segments = segmentation(gnome,t1);
    // Loop this if you want recursive.  Otherwise single-pass.
    std::cout << segments.size() << " Segments Generated Across " << gnome.length() << " Nucleotides" << std::endl;
    std::cout << "Stage 1 Clustering..." << std::endl;
    std::vector<cluster> stage1segments = clusterzontal(segments,0);
    int stagesize = 5;
    std::vector<cluster> stage2segments;
    std::vector<cluster> oristage2segments;
    while(stagesize != stage2segments.size()) {
        stagesize = stage2segments.size();
        stage2segments = clusterzontal2(stage1segments,t2);
        stage1segments = stage2segments;
    }
    // Also single-pass.
    std::cout << stage2segments.size() << " Clusters Remaining Post-Merge" << std::endl;
    oristage2segments = stage2segments;
    for(int i=0;i<range.size();i++) {
        std::unordered_set<std::string> dmjsdlib;
        stage2segments = oristage2segments;
        std::cout << "Running T3 = " << range[i] << std::endl;
        outfile = argv[4];
        outfile = outfile + "." + range[i];
        t3 = std::stod(range[i]);
        std::vector<cluster> stage3segments;
        std::cout << "Stage 2 Clustering..." << std::endl;
        stage3segments = clusterbeyondX(stage2segments,t3);
        std::cout << stage3segments.size() << " Clusters Remaining Post-Merge" << std::endl;
        std::cout << "Printing Output to " << outfile << std::endl;
        fasfilePrinter(outfile,stage3segments);
    }
    return 0;
    
}
