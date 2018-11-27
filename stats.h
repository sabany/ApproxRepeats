/*
 *  Methods for assessing the statistical signifiance of molecular
 *  sequence features by using general scoring schemes, writting by
 *  Samuel Karlin and F. Atschul in 1989-90.
 */

#ifndef __STATS_H_INCLUDED__
#define __STATS_H_INCLUDED__

#define ABS(a)              ((a)>(0)?(a):(-(a)))
#define MEGA                (1024*1024)
#define SENSIBILITY_LAMBDA  1e-9
#define DBL_MAX             10e16
#define NBLOOPS_K           20
#define DELTA_K             1

using namespace std;

int32_t substitutionMatrix [4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
intType NbLetters          [2][4] = {{0,0,0,0},{0,0,0,0}};   // number of 'A','T','G','C' for each file
double ** freqLetters = NULL;                                // 2x4, frequency of 'A','T','G','C' for each file
double ** freqBackground = NULL;                             // 4x4

int32_t MinScore(double K, double Lambda, intType ref_size, intType query_size, double MaxEvalue)
{
    return (int32_t) floor((1 / Lambda) * (log(K * ref_size * query_size) - log(MaxEvalue)));
}

double ** dbl_directtable(uint32_t i, uint32_t j)
{
    double ** itable;
    double *  jtable;
    uint32_t k;
    
    itable = new double*[i + (i * j)];
    jtable = new double[i * j];

    jtable = (double *) (itable + i);
    
    for(k=0 ; k<i; k++)
        itable[k] = jtable + k * j;
    
    return itable;
}

//  Compute letters frequency
double ** computeLettersFrequency(intType nb_letters[2][4])
{
    intType nbL[2];
    double ** letters_frequency;
    
    nbL[0] =  nb_letters[0][0] + nb_letters[0][1] + nb_letters[0][2] + nb_letters[0][3];
    nbL[1] =  nb_letters[1][0] + nb_letters[1][1] + nb_letters[1][2] + nb_letters[1][3];
    
    if (nbL[0] == 0) {
        nb_letters[0][0] = nb_letters[0][1] = nb_letters[0][2] = nb_letters[0][3] = 1;
        nbL[0] = 4;
    }
    
    if (nbL[1] == 0) {
        nb_letters[1][0] = nb_letters[1][1] = nb_letters[1][2] = nb_letters[1][3] = 1;
        nbL[1] = 4;
    }
    
    /* table allocation */
    letters_frequency = dbl_directtable(2,4);
    
    /* freq computed */
    {
        uint32_t i,j;
        for(i=0;i<2;i++)
            for(j=0;j<4;j++)
                letters_frequency[i][j] = (double) nb_letters[i][j] / nbL[i];
    }
    return letters_frequency;
}

// Compute mutation probabilty on single nucleotides
double ** computeBackgroundFrequency(double ** letters_frequency /* [2][4] */)
{
    
    double ** background_frequency;
    
    /* table allocation */
    background_frequency = dbl_directtable(4,4);
    {
        uint32_t i,j;
        for(i=0;i<4;i++)
            for(j=0;j<4;j++)
                background_frequency[i][j] =  letters_frequency[0][i] * letters_frequency[1][j];
    }
    
    /* 3) return parametre */
    return background_frequency;
}

// Compute Lambda
#define COMPUTELAMBDA(lambda)\
{\
    uint32_t i,j;\
    S = 0.0;\
    for (i=0; i<4; i++)\
        for (j=0; j<4; j++)\
            S += (freq_background[i][j]) * exp(lambda*(double)(substitutionMatrix[i][j]));\
}

double computeLambda(double ** freq_background /* [4][4] */ ) {
    
    double lambda_lower = 0.0;
    double lambda_upper = 1.0;
    double lambda = 0.0;
    double S = 0.0;
    
    /*1) check feasability  */
    {
        uint32_t i,j;
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                S += (freq_background[i][j]) * (double) (substitutionMatrix[i][j]);
    }
    
    if (S >= 0.0) {
        cout << "ERROR: Lambda value cannot be computed:" << endl;
        cout << "\tA common reason is a strong AT/GC sequence bias, so please fix it with" << endl;
        cout << "\tthe \'-Lambda\' option, or change the scoring system with the \'-s\' option." << endl;
        return 0;
    }
    
    /*2) looking for lambda upper */
    COMPUTELAMBDA(lambda_upper);
    while (S < 1) {
        lambda_lower = lambda_upper;
        lambda_upper *= 2.0;
        COMPUTELAMBDA(lambda_upper);
    }
    
    
    /*3) looking for lambda exactly */
    while (lambda_upper - lambda_lower > SENSIBILITY_LAMBDA) {
        lambda =  ( lambda_upper + lambda_lower ) * 0.5;
        COMPUTELAMBDA(lambda);
        if (S > 1.0) {
            lambda_upper = lambda;
        } else {
            lambda_lower = lambda;
        }
    }
    return lambda;
}

// compute K
double computeK(double ** freq_background /* [4][4] */, double lambda) {
    
    double *pr = NULL, *pr_new = NULL, *pr_tmp = NULL;
    double denominator = 0.0;
    double numerator = 0.0;
    double Pb = 0.0 , Ek = 0.0 , C = 0.0;
    int32_t i=0, k=0, scoreMin = 0, scoreMax = 0;
    int32_t lengthProb, zeroPosition;
    
    
    /*
     * Tests values from BLAST:
     */
    
    /*
     * Blast values that must be checked
     *  -C +1/-3 (and 50 %GC)
     *  lambda = 1.371;
     *  K      = 0.711;
     */
    
    /* Find minScore,MaxScore to compute the tab length */
    {
        uint32_t i,j;
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                if ( substitutionMatrix[i][j] > scoreMax )
                    scoreMax = substitutionMatrix[i][j];
                if ( substitutionMatrix[i][j] < scoreMin)
                    scoreMin =  substitutionMatrix[i][j];
            }
        }
    }
    
    lengthProb   = (scoreMax - scoreMin + 1) * NBLOOPS_K + 1;
    zeroPosition = ABS(scoreMin) * NBLOOPS_K;
    
    
    /*1) initialise vector probability tables */
    if (lengthProb > 4 * MEGA) {
        cout<<"ERROR: compute K is not possible, because the scores given (-s option) are too large"<<endl;
        return 0;
    }
    
    try{
        pr = new double[lengthProb];
    }catch(std::bad_alloc& ba){
        std::cerr << "computeK _ bad_alloc caught: " << ba.what() <<endl;
    }
    try{
        pr_new = new double[lengthProb];
    }catch(std::bad_alloc& ba){
        std::cerr << "computeK _ bad_alloc caught: " << ba.what() <<endl;
    }
    
    for (i = 0; i < lengthProb; i++) {
        pr[i] = 0.0;
        pr_new[i] = 0.0;
    }
    
    /* Initial step for K = 0 */
    pr[zeroPosition] = 1.0;
    
    /*
     * 2) Compute denominator lambda*E[ S_1 * e^(lambda*S_1)]
     */
    denominator = 0.0;
    {
        uint32_t i,j;
        for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                denominator += freq_background[i][j] * substitutionMatrix[i][j] * exp(lambda * substitutionMatrix[i][j]);
            }
        }
    }
    denominator *= lambda;
    
    
    /*
     * 3) numerator 1/k*SUM (E[e^lambda*S_k;S_k<0]+Prob(S_k>=0),k,1,SENSIBILITY_K)
     */
    
    for (k = 1; k < NBLOOPS_K; k++) {
        int32_t x;
        
        /*
         * 3.1) compute a new probability
         */
        
        for (x = ABS(scoreMin); x < lengthProb - ABS(scoreMax); x++){
            uint32_t i,j;
            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++) {
                    pr_new[x +  substitutionMatrix[i][j]] += pr[x] * freq_background[i][j];
                }
            }
        }
        
        
        /*
         * 3.2) permutation pr and pr_new
         */
        pr_tmp = pr;
        pr = pr_new;
        pr_new = pr_tmp;
        
        for (i = 0; i < lengthProb; i++)
            pr_new[i] = 0.0;
        
        
        Pb = 0.0;		/* ProbabilitÈ Prob(S_k) */
        Ek = 0.0;		/* Esperance   E[e^(lambda*S_k);S_k<0 */
        
        
        /*
         * 3.4) esperance on the negative side (<0)
         */
        for (i = 0; i < zeroPosition; i++) {
            int32_t Sk = i - zeroPosition;
            
            if (exp(lambda * (double)Sk) > DBL_MAX) {
                cout << "ERROR: compute K,  exponent is too large" << endl;
                return 0;
            }
            Ek += pr[i] * exp(lambda * (double)Sk);
        }
        
        /*
         * 3.5) probabilities on the postif sidee (>=0)
         */
        for (i = zeroPosition ; i < lengthProb; i++)
            Pb += pr[i];
        
        numerator += (Pb + Ek) / (double) k;
    }
    
    C = exp(-2 * numerator) / (denominator);
    
    free(pr);
    free(pr_new);
    
    return C * (lambda * DELTA_K) / (1 - exp(-lambda * DELTA_K));
}

void fillSubstitutionMatrix() {
    uint32_t i, j;
    
    for(i=0; i<4 ; i++)
        for(j=0; j<4; j++)
            if (i==j)
                substitutionMatrix[i][j] = commonData::matchScore;
            else
                substitutionMatrix[i][j] = commonData::misMatScore;
}

void fillNbLetters(seqFileReadInfo & seqFile, uint32_t refORquery) {  // refORquery=0 denotes Reference, refORquery=1 denotes Query
    NbLetters[refORquery][0] = seqFile.numOfAs;;
    NbLetters[refORquery][1] = seqFile.numOfCs;
    NbLetters[refORquery][2] = seqFile.numOfGs;
    NbLetters[refORquery][3] = seqFile.numOfTs;
}

bool setMinScore(intType totalRBases, intType totalQBases) {
    
    freqLetters = computeLettersFrequency(NbLetters);
    freqBackground = computeBackgroundFrequency(freqLetters);
    
    commonData::lambdaBlast = computeLambda(freqBackground);
    if (!commonData::lambdaBlast)
        return false;
    
    commonData::kBlast = computeK(freqBackground,commonData::lambdaBlast);
    if (!commonData::kBlast)
        return false;
    
    commonData::minScore = MinScore(commonData::kBlast, commonData::lambdaBlast, totalRBases, totalQBases, commonData::expectationValue);
    return true;
}

#endif
