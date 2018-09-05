/*
    Title:            MrModeltest2
    Version:          UNIX, MacOSX, Win32
    Latest changes:   Tue 04 Sep 2018
    Programmer:       Johan Nylander
                      Uppsala University
                      E-mail: johan.nylander@ebc.uu.se
    Notes:            This is a simplified version of David Posada's modeltest v.3.6.
                      The difference between MrM. and M. is that MrM. only tests 24
                      substittion models; the 24 nucleotide-substitution models
                      currently implemented in MrBayes v3. Modeltest tests 56 models,
                      some of them not implemeted in MrBayes.
                      In addition, MrM allow the user to choose different hierarchies
                      for the likelihood-ratio tests, and prints an example of how to
                      implement the selected model in MrBayes v3.
    Credits:          David Posada are truly thanked for supplying the excellent code in Modeltest!
    Version history:  MrModeltest2.2, 2005-02-01: Fixed a bug affecting
                      the printing of PAUP and MrBayes blocks.
                      MrModeltest2.3 2008-04-16: The MrBayes block was not printed correctly:
                      The settings made by the Prset command was overwritten by the Lset
                      command if issued before the Lset command. Thanks to Ted Schultz.
                      2016-11-02: Took care of some compiler warnings and cleaned the code.
                      MrModeltest2.4 2018-09: The output from PAUP* changed. MrModeltest2.4
                      can now only read the new format. Version 2.4 does not report model
                      averaging output. Thanks to Leila Carmona.
    ===========================================================================
            Below are David's original notes on Modeltest.
    ===========================================================================
    Title: modeltest
    Programmer: David Posada
    Date started: 03/05/98
    Purpose: Evaluate the fit of several models of evolution to a given data and unrooted tree.
    Compare different models of DNA substitution using likelihood ratio test (and a chi-square distribution)
    or/and the AIC criterion

    COPYRIGHT
    --------
    Copyright (C) 1998-2004 David Posada
    Facultad de Biolog’a, Universidad de Vigo. 36200 Vigo, Spain
    dposada@uvigo.es

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

    HISTORY
    -------
    Version 2.0 (June 99):     Models TrN, K81 and its equivalents with equal base frequencies added
    Version 2.1 (October 99):  Changed (-1) #free parameters (eg. JC has 0 instead of 1).
    Version 3.0beta1 (Nov 99): Models TIM, TIMef, TVM and TVMef added. Now we have 56 models.
    Version 3.0beta3 (Jan 00): Now PAUP* outputs base frequencies estimates in the likelihood scores file,
                               so changes were done to read this new file
    Version 3.0 (February 00): Several aesthetic changes
    Version 3.01(March 00):    The program ouputs now a block of PAUP commands
                               to implement the likelihood settings for the best-fit model
    Version 3.02 (June 00):    The frequencies of C and G were interchanged in the output. Fixed. Thanks to Carlos Lopez-Vaamonde
    Version 3.03 (June 00):    The mixed chi-square distribution is added for I and G tests
    Version 3.04 (July 00):    Several cosmetic changes
    Version 3.05 (Feb 01):     In the windows version, the AIC[55] gave an AIC of 0 to the GTRIG. Now AIC[56]. (Juan Suarez)
                               TIM+G reported invariable sites instead of gamma shape (Cymon Cox)
    Version 3.06 (Apr 01):     Print likelihood scores by default
                               In the windows version there was a bug by which the file scores.txt was always
                               the standard input (Andy Vierstraete)
                               Using GNU licencese (I should have done this a long time ago) (thanks to Naoki Takebayashi)
    Version 3.07 (Apr 01):     Mistake in the print likelihoods corrected
                               Some model names were not displayed complete in the screen (e.g TrNef+I+G)
    Version 3.1 (March 02):    Calculate delta AIC and Akaike weights
                               Removed option to turn off the use of mixed chi square for I and G LRTs
    Version 3.2 (March 03):    Aesthetic changes
                               TrN+I had 5 df instead of 6
    Version 3.3 (Nov 03):      Option to include branch length estimates as parameters
                               Option to calculate AICc
                               Changing some option letters accordingly
    Version 3.4 (March 04):    There was a typo printing the Rd value for K81uf+I.
                               It was printing p-inv instead. (thanks to Michael Sorenson)
    Version 3.5 (May 04):      This is a minor update that does not affect the calculations. AIC weights were sorted by
                               their value, but because these can be almost zero (zero for the computer) for several models,
                               their order would not make sense in the light of the AIC values.
                               Now the program order the AIC weights by the AIC scores.
    Version 3.6 beta (Nov 04): The program includes now model averaged estimates
                               The program includes now variable importance calculations
                               New option (-w) to define confidence interval of models used for model-averaged estimates. By
                               default this interval is 1.0, so all models are included in model-averaged estimates
                               Use double precision in the AIC calculator (thanks to Renee Park)
                               Likelihoods and number of parameters are now printed in the AIC table
                               Aesthetic changes
                               Argument for specifying sample size is now -n (it was -c)
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>

/* Constants */
#define BIGX           20.0                           /* max value to represent exp (x) */
#define LOG_SQRT_PI    0.5723649429247000870717135    /* log (sqrt (pi)) */
#define I_SQRT_PI      0.5641895835477562869480795    /* 1 / sqrt (pi) */
#define Z_MAX          6.0                            /* maximum meaningful z value */
#define ex(x)          (((x) < -BIGX) ? 0.0 : exp (x))
#define MAX_PROB       0.999999
#define MIN_PROB       0.000001
#define PROGRAM_NAME   "MrModeltest"
#define VERSION_NUMBER "2.4"
#define SUCCESS        1
#define FAILURE        0
#define YES            1
#define NO             0
#define BIGNUMBER      9999999
#define NA             -99999
#define NUM_MODELS     24

/* Structures */
typedef struct {
    float ln;
    int parameters;
    char *name;
} ModelSt;

/* Prototypes */
static void ReadArgs(int, char**);
static void RecognizeInputFormat();
static void ReadPaupScores();
static void Initialize();
static void ReadScores();
static void PrintTitle(FILE *fp);
static void PrintDate(FILE *fp);
static void RatioCalc();
static void hLRT();
static void hLRT2();
static void hLRT3();
static void hLRT4();
static void CalculateAIC();
static void AkaikeWeights ();
static int AICCalc();
static void AICfile();
static void SetModel(char[]);
static void Allocate();
static void Free();
static void PrintUsage();
static void HLRTAttention(char *first, char *second, char *third, char *fourth);
static void Output(char *selection, float value);
static void PrintPaupBlock(int ishLRT);
static void PrintMbBlock(int ishLRT);
static double LRT(ModelSt *model0, ModelSt *model1);
static double LRTmix(ModelSt *model0, ModelSt *model1);
static void PrintRunSettings();
/*static void ModelAveraging();*/
/*static void AverageEstimates (char *parameter, int numModels, int modelIndex[], int estimateIndex[],
double *importance, double *averagedEstimate, double minWeightToAverage);*/
/*static double FindMinWeightToAverage ();*/
/* static char *CheckNA (double value);*/
float ChiSquare (float x, int);
float Normalz (float);
float TestEqualBaseFrequencies(ModelSt *, ModelSt *);
float TestTiequalsTv(ModelSt *, ModelSt *);
float TestEqualTiAndEqualTvRates(ModelSt *, ModelSt *);
float TestEqualSiteRates(ModelSt *, ModelSt *);
float TestInvariableSites(ModelSt *, ModelSt *);

/* Global variables */
ModelSt *JC, *F81;
ModelSt *JC, *JCI, *JCG, *JCIG, *F81, *F81I, *F81G, *F81IG, *K80, *K80I, *K80G, *K80IG;
ModelSt *HKY, *HKYI, *HKYG, *HKYIG, *SYM, *SYMI, *SYMG, *SYMIG, *GTR, *GTRI, *GTRG, *GTRIG;
float score[176]; /* Check this length. 175? Previous:: 168 */
float ln[NUM_MODELS];
float AIC[NUM_MODELS];
float wAIC[NUM_MODELS];
int orderedAIC[NUM_MODELS];
float alpha, minAIC;
int format, file_id;
int print_scores;
int DEBUGLEVEL;
char *modelhLRT;
char *modelhLRT2;
char *modelhLRT3;
char *modelhLRT4;
char *modelAIC;
char *nexttok1;
char *nexttok2;
char filename[80];
char infilename[80];
ModelSt *model;
ModelSt *order;
FILE *fpin;
int mixchi;
int numTaxa, numBL, sampleSize;
int useAICc, useBL;
int usehLRT2 = NO;
int usehLRT3 = NO;
int usehLRT4 = NO;
int lastModelConfidence;
/*float averagingConfidenceInterval;*/
double cumConfidenceWeight;

/* Parameter estimates for the selected model */
float fA, fC, fG, fT;
float TiTv;
/*float Ra, Rb, Rc, Rd, Re, Rf;*/
float rAC, rAG, rAT, rCG, rCT, rGT;
float shape;
float pinv;
float theln;

/****************************** MAIN ***********************************/
int main(int argc, char **argv)
{
    float start, secs;
    Allocate();
    alpha = 0.01;                      /* default level of significance (aprox Bonferroni)  */
    mixchi = YES;                      /* by default use mixed chi-square distribution */
    print_scores = YES;                /* by default print the likelihood scores for all models */
    numTaxa = 0;
    sampleSize = 0;
    useBL = NO;                        /* by default do not include branch length estimates as parameters */
    useAICc = NO;                      /* by default do not use the AICc correction */
    /*averagingConfidenceInterval = 1.0;*/ /* by default include all models in model-averaged estimates */

    start = clock();
    ReadArgs(argc, argv);
    file_id = isatty(fileno(stdin));
    if (file_id) {
        fprintf(stderr, "\n\nNo input file\n\n");
        PrintUsage();
        /*getchar();*/    /*For windows.*/
        exit(1);
    }
    PrintTitle(stdout);
    PrintDate(stdout);
    RecognizeInputFormat();
    PrintRunSettings();

    /* Do hLRTs */
    printf("\n\n\n\n---------------------------------------------------------------");
    printf("\n*                                                             *");
    printf("\n*         HIERARCHICAL LIKELIHOOD RATIO TESTS (hLRTs)         *");
    printf("\n*                                                             *");
    printf("\n---------------------------------------------------------------\n");

    if (usehLRT4 == YES) {
        hLRT4();
        modelhLRT = modelhLRT4;
    }
    else if (usehLRT3 == YES) {
        hLRT3();
        modelhLRT = modelhLRT3;
    }
    else if (usehLRT2 == YES) {
        hLRT2();
        modelhLRT = modelhLRT2;
    }
    else {
        hLRT();
        printf("\n\n\n ** Hierarchical Likelihood Ratio Tests (using hLRT2) **\n");
        hLRT2();
        printf("\n\n\n ** Hierarchical Likelihood Ratio Tests (using hLRT3) **\n");
        hLRT3();
        printf("\n\n\n ** Hierarchical Likelihood Ratio Tests (using hLRT4) **\n");
        hLRT4();
    }

    SetModel(modelhLRT);

    if (format == 0) {
        if (shape > 999) { /* alpha shape = infinity */
            printf("\n\nWARNING: Although the model %s was initially selected, gamma (G) was removed ",modelhLRT);
            printf("because the estimated shape equals infinity, which implies equal rates among sites.");
            /* removing +G */
            nexttok1 = strtok(modelhLRT, "+");
            strcpy(modelhLRT, nexttok1);
            if (!strcmp(strtok(NULL, "+"),"I")) {
                strcat(modelhLRT,"+I");
            }
        }
        else if (usehLRT2 == usehLRT3 && usehLRT2 == usehLRT4 && usehLRT3 == usehLRT4) { /* Old: if (usehLRT2 == usehLRT3 == usehLRT4) */
            Output (modelhLRT,0);
            PrintPaupBlock (YES);
            PrintMbBlock (YES);
        }
        else {
            HLRTAttention (modelhLRT, modelhLRT2, modelhLRT3, modelhLRT4);
            Output (modelhLRT,0);
            PrintPaupBlock (YES);
            PrintMbBlock (YES);
        }
    }
    else {
        printf("\n hLRT model = %s", modelhLRT);
        if (modelhLRT2) {
            printf("\n hLRT2 model = %s", modelhLRT2);
        }
        if (modelhLRT3) {
            printf("\n hLRT3 model = %s", modelhLRT3);
        }
        if (modelhLRT4) {
            printf("\n hLRT4 model = %s", modelhLRT4);
        }
    }

    /* Do AIC */
    if (useAICc == YES) {
        printf("\n\n\n\n\n---------------------------------------------------------------");
        printf("\n*                                                             *");
        printf("\n*        SECOND ORDER AKAIKE INFORMATION CRITERION (AICc)        *");
        printf("\n*                                                             *");
        printf("\n---------------------------------------------------------------\n");
    }
    else {
        printf("\n\n\n\n\n---------------------------------------------------------------");
        printf("\n*                                                             *");
        printf("\n*             AKAIKE INFORMATION CRITERION (AIC)              *");
        printf("\n*                                                             *");
        printf("\n---------------------------------------------------------------\n");
    }
    CalculateAIC();
    SetModel(modelAIC);
    if (format == 0) {
        Output (modelAIC,minAIC);
        PrintPaupBlock(NO);
        PrintMbBlock(NO);
    }
    else if (useAICc == YES) {
        printf("\n AICc model = %s", modelAIC);
    }
    else {
        printf("\n AIC model = %s", modelAIC);
    }
    AkaikeWeights();
    /*ModelAveraging();*/
    secs = (double)(clock() - start) / CLOCKS_PER_SEC;
    Free();
    printf("\n\n_________________________________________________________________________");
    printf("\nTime processing: %G seconds", secs);
    printf("\nIf you need help type '-?' or '-h' in the command line of the program");
    fprintf(stderr, "\nProgram is done.\n\n");
    return 0;
}

/******************** PrintRunSettings **************************/
static void PrintRunSettings()
{
    int i;
    /* Check settings */
    fprintf (stdout,"\n\nRun settings\n");
    if (sampleSize > 0) {
        useAICc = YES;
        /* Check the data is large enough*/
        if (sampleSize <= model[NUM_MODELS-1].parameters) {
            fprintf (stdout, "\n\nYou have more parameters than data for some models!");
            fprintf (stdout, "\nCalculations cannot be performed. Exiting the program ...\n");
            exit(1);
        }
        fprintf (stdout, "\n Using the AICc correction");
        fprintf (stdout, "\n   sample size = %d", sampleSize);
    }
    else {
        useAICc = NO;
        fprintf (stdout, "\n Using the standard AIC (not the AICc)");
    }
    if (numTaxa > 0) {
        useBL = YES;
        fprintf (stdout, "\n Using branch lengths as parameters");
        fprintf (stdout, "\n   number of taxa = %d (%d branch lengths)", numTaxa, numBL);
        for (i = 0; i < NUM_MODELS; i++) {
            model[i].parameters += numBL;
            order[i].parameters += numBL;
        }
    }
    else {
        useBL = NO;
        fprintf (stdout, "\n Not using branch lengths as parameters");
    }
    if (usehLRT4 == YES) {
        fprintf (stdout, "\n Printing results based on the hLRT4 hierarchy");
    }
    else if (usehLRT3 == YES) {
        fprintf (stdout, "\n Printing results based on the hLRT3 hierarchy");
    }
    else if (usehLRT2 == YES) {
        fprintf (stdout, "\n Printing results based on the hLRT2 hierarchy");
    }
    else {
        fprintf (stdout, "\n Running all four hierarchies for the hLRT");
        fprintf (stdout, "\n Printed parameter values are from the hLRT1 hierarchy");
    }
}

/******************** ReadArgs **************************/
static void ReadArgs(int argc,char **argv)
{
    int i;
    char flag;

    for (i = 1; i < argc; i++) {
        argv[i]++;
        flag = *argv[i];
        argv[i]++;
        switch (flag) {
        case 'd':
            DEBUGLEVEL = atoi(argv[i]);
            if (DEBUGLEVEL >= 2)
                fprintf(stderr,"DEBUGLEVEL set to %d\n", DEBUGLEVEL);
            break;
        case 'a':
            alpha = atof(argv[i]);
            fprintf(stderr, "alpha level set to %.6f\n", alpha);
            break;
        case 'n':
            sampleSize = atoi(argv[i]);
            break;
        case 't':
            numTaxa = atoi(argv[i]);
            numBL = 2 * numTaxa - 3;
            break;
        case 'l':
            printf("\n LRT CALCULATOR MODE \n");
            RatioCalc();
            break;
        case 'i':
            printf("\n AIC CALCULATOR MODE \n");
            AICCalc();
            break;
        case 'f':
            AICfile();
            break;
            case '?':
            PrintUsage();
            exit(1);
            break;
        case 'h':
            PrintUsage();
            exit(1);
            break;
        case '2':
            usehLRT2 = YES;
            break;
        case '3':
            usehLRT3 = YES;
            break;
        case '4':
            usehLRT4 = YES;
            break;
        case 'v':
            fprintf(stderr,"%s version %s\n",PROGRAM_NAME, VERSION_NUMBER);
            exit(1);
            break;
//        case 'w':
//            averagingConfidenceInterval = atof(argv[i]);
//            if (averagingConfidenceInterval > 1) {
//                fprintf (stderr, "\nError: confidence interval cannot be > 1");
//                exit (1);
//            }
//            else if (averagingConfidenceInterval <= 0) {
//                fprintf (stderr, "\nError: confidence interval cannot be <= 0");
//                exit (1);
//            }
//            break;
        default:
            fprintf(stderr,"Unknown argument on the command line '%c'\n",flag);
            exit(1);
            break;
        }
    }
}

/********************* RecognizeInputFormat ***********************/
static void RecognizeInputFormat()
{
    int iochar;

    if (stdin == NULL) {
        fprintf(stderr,"Error opening the input file");
        exit(0);
    }
    iochar = getc(stdin);
    if (iochar == (int)'T') {   /* In the Paup matrix, in the first line there is the word 'Tree'*/ /*Att göra: kolla input format? */
        ungetc(iochar,stdin);
        printf("\nInput format: Paup matrix file \n");
        format = 0;
        ReadPaupScores();
    }
    else {
        ungetc(iochar,stdin);
        printf("\nInput format: raw log likelihood scores \n");
        format = 1;
        ReadScores();
    }
}

/***************************** ReadPaupScores ********************************/
static void ReadPaupScores()
{
    int iochar;
    int i,j,k;
    char string [120];
    i=0;

    while (!feof(stdin)) {
        iochar = getc(stdin);
        if (isdigit(iochar)) {
            ungetc(iochar,stdin);
            scanf("%f", &score[i]);
            if (DEBUGLEVEL >= 2) {
                fprintf(stdout,"\nINFO:   Storing %f in score[%d]", score[i],i);
            }
            i++;
        }
        if (isalpha (iochar)) {
            ungetc(iochar,stdin);
            scanf("%s",string);
            if (DEBUGLEVEL >= 2) {
                fprintf(stdout,"\nINFO:   Reading string %s", string);
            }
            if (strcmp(string,"infinity") == 0) {
                score[i] = 999.999;
                if (DEBUGLEVEL >= 2)
                fprintf(stdout,"INFO:   Storing %f in score[%d]", score[i],i);
                i++;
            }
        }
    }
    if (ferror(stdin)) {
        perror ("MrModeltest");
        clearerr(stdin);
    }
    Initialize();
    if(print_scores == YES) {
        printf("\n\n** Log Likelihood scores **");
        printf("\n%-12.12s\t\t\t+I\t\t+G\t\t+I+G", " ");
        for(k = 0; k < NUM_MODELS; k += 4) {
            printf("\n%-10.10s =\t%9.4f\t%9.4f\t%9.4f\t%9.4f", model[k].name, model[k].ln, model[k+1].ln, model[k+2].ln, model[k+3].ln);
        }
        printf("\n\n");
    }
    for (j = 0; j < NUM_MODELS; j++) {
        if (model[j].ln == 0 || i < 175) { /*Att göra: kolla detta värde! Ev 175. */
            printf("\n\nThe input file is incomplete or incorrect. \nAre you using the most updated block of PAUP* commands?. ");
            printf("\nThis version of MrModeltest is not compatible with versions of PAUP* older than PAUP*4.0beta3");
            printf("\nCheck the Modeltest and PAUP* web pages");
            exit(0);
        }
     }
}

/************** Initialize. Modified by Johan 2002-03-18 **********************/
/* New versions of paup (at least year 2018) prints a different output.*/
/* Last edits Tue 04 Sep 2018 02:42:12 PM CEST */
void Initialize()
{
    JC =     model;      model[0].ln = order[0].ln = score[1]; /*1*/
    JCI =    model + 1;  model[1].ln = order[1].ln = score[3]; /*3*/
    JCG =    model + 2;  model[2].ln = order[2].ln = score[6]; /*6*/
    JCIG =   model + 3;  model[3].ln = order[3].ln = score[9]; /*9*/
    F81 =    model + 4;  model[4].ln = order[4].ln = score[13]; /*13*/
    F81I =   model + 5;  model[5].ln = order[5].ln = score[19]; /*19*/
    F81G =   model + 6;  model[6].ln = order[6].ln = score[26]; /*26*/
    F81IG =  model + 7;  model[7].ln = order[7].ln = score[33]; /*33*/
    K80 =    model + 8;  model[8].ln = order[8].ln = score[41]; /*41*/
    K80I =   model + 9;  model[9].ln = order[9].ln = score[44]; /*44*/
    K80G =   model + 10; model[10].ln = order[10].ln = score[48]; /*48*/
    K80IG =  model + 11; model[11].ln = order[11].ln = score[52]; /*52*/
    HKY =    model + 12; model[12].ln = order[12].ln = score[57]; /*57*/
    HKYI =   model + 13; model[13].ln = order[13].ln = score[64]; /*64*/
    HKYG =   model + 14; model[14].ln = order[14].ln = score[72]; /*72*/
    HKYIG =  model + 15; model[15].ln = order[15].ln = score[80]; /*80*/
    SYM =    model + 16; model[16].ln = order[16].ln = score[89]; /*89*/
    SYMI =   model + 17; model[17].ln = order[17].ln = score[97]; /*96*/
    SYMG =   model + 18; model[18].ln = order[18].ln = score[106]; /*104*/
    SYMIG =  model + 19; model[19].ln = order[19].ln = score[115]; /*112*/
    GTR =    model + 20; model[20].ln = order[20].ln = score[125]; /*121*/
    GTRI =   model + 21; model[21].ln = order[21].ln = score[137]; /*132*/
    GTRG =   model + 22; model[22].ln = order[22].ln = score[150]; /*144*/
    GTRIG =  model + 23; model[23].ln = order[23].ln = score[163]; /*156*/
    /* free parameters */
    /* JC */
    model[0].parameters = order[0].parameters = 0;
    model[1].parameters = order[1].parameters = 1;
    model[2].parameters = order[2].parameters = 1;
    model[3].parameters = order[3].parameters = 2;
    /* F81*/
    model[4].parameters = order[4].parameters = 3;
    model[5].parameters = order[5].parameters = 4;
    model[6].parameters = order[6].parameters = 4;
    model[7].parameters = order[7].parameters = 5;
    /* K80 */
    model[8].parameters = order[8].parameters = 1;
    model[9].parameters = order[9].parameters = 2;
    model[10].parameters = order[10].parameters = 2;
    model[11].parameters = order[11].parameters = 3;
    /* HKY */
    model[12].parameters = order[12].parameters = 4;
    model[13].parameters = order[13].parameters = 5;
    model[14].parameters = order[14].parameters = 5;
    model[15].parameters = order[15].parameters = 6;
    /* SYM */
    model[16].parameters = order[16].parameters = 5;
    model[17].parameters = order[17].parameters = 6;
    model[18].parameters = order[18].parameters = 6;
    model[19].parameters = order[19].parameters = 7;
    /* GTR */
    model[20].parameters = order[20].parameters = 8;
    model[21].parameters = order[21].parameters = 9;
    model[22].parameters = order[22].parameters = 9;
    model[23].parameters = order[23].parameters = 10;
    /* */
    order[0].name =  "JC";
    order[1].name =  "JC+I";
    order[2].name =  "JC+G";
    order[3].name =  "JC+I+G";
    order[4].name =  "F81";
    order[5].name =  "F81+I";
    order[6].name =  "F81+G";
    order[7].name =  "F81+I+G";
    order[8].name =  "K80";
    order[9].name =  "K80+I";
    order[10].name = "K80+G";
    order[11].name = "K80+I+G";
    order[12].name = "HKY";
    order[13].name = "HKY+I";
    order[14].name = "HKY+G";
    order[15].name = "HKY+I+G";
    order[16].name = "SYM";
    order[17].name = "SYM+I";
    order[18].name = "SYM+G";
    order[19].name = "SYM+I+G";
    order[20].name = "GTR";
    order[21].name = "GTR+I";
    order[22].name = "GTR+G";
    order[23].name = "GTR+I+G";
    /* */
    model[0].name =  "JC";
    model[1].name =  "JC+I";
    model[2].name =  "JC+G";
    model[3].name =  "JC+I+G";
    model[4].name =  "F81";
    model[5].name =  "F81+I";
    model[6].name =  "F81+G";
    model[7].name =  "F81+I+G";
    model[8].name =  "K80";
    model[9].name =  "K80+I";
    model[10].name = "K80+G";
    model[11].name = "K80+I+G";
    model[12].name = "HKY";
    model[13].name = "HKY+I";
    model[14].name = "HKY+G";
    model[15].name = "HKY+I+G";
    model[16].name = "SYM";
    model[17].name = "SYM+I";
    model[18].name = "SYM+G";
    model[19].name = "SYM+I+G";
    model[20].name = "GTR";
    model[21].name = "GTR+I";
    model[22].name = "GTR+G";
    model[23].name = "GTR+I+G";
}

/******************* ReadScores ************************/
static void ReadScores()
{
    float score[NUM_MODELS];
    int i,j;
    i=0;

    score[NUM_MODELS-1] = 0;
    while (!feof(stdin)) {
        scanf("%f", &score[i]);
        i++;
    }
    Initialize();
    for (j = 0; j < NUM_MODELS; j++) {
        if (model[j].ln == 0 || i < 175) {
            printf("\n\nThe input file is incomplete or incorrect. \nAre you using the most updated block of PAUP* commands?.\n ");
            exit(0);
        }
    }
    if (i > NUM_MODELS + 1) {
        printf("\n\nThe input file has more than %d scores", NUM_MODELS);
        exit(0);
    }
}

/******************* LRT ******************************/
/* performs a likelihood ratio test */
static double LRT(ModelSt *model0, ModelSt *model1)
{
    double delta;
    double prob;
    int df;

    delta = 2 * (model0->ln - model1->ln);
    df = model1->parameters - model0->parameters;
    if (delta == 0) {
        prob = 1.0;
    }
    else {
        prob = ChiSquare(delta, df);
    }
    printf("\n   Null model = %-9.9s\t\t  -lnL0 = %.4f", model0->name, model0->ln);
    printf("\n   Alternative model = %-9.9s\t  -lnL1 = %.4f", model1->name, model1->ln);
    printf("\n   2(lnL1-lnL0) = %9.4f\t\t      df = %d ", delta, df);
    if (prob == 1.0) {
        printf("\n   P-value = >%f", MAX_PROB);
    }
    else if (prob < 0.000001) {
        printf("\n   P-value = <%f", MIN_PROB);
    }
    else {
        printf("\n   P-value =  %f", prob);
    }
    return prob;
}

/******************* LRTmix ******************************/
/* performs a likelihood ratio test and uses mixed chi2*/
static double LRTmix(ModelSt *model0, ModelSt *model1)
{
    double delta, prob;
    int df;

    delta = 2 * (model0->ln - model1->ln);
    df = model1->parameters - model0->parameters;
    if (delta == 0) {
        prob = 1.0;
    }
    else {
        if (df == 1) {
            prob = ChiSquare(delta,df)/2;
        }
        else {
            prob = (ChiSquare(delta,df-1) + ChiSquare(delta,df)) / 2;
        }
    }
    printf("\n   Null model = %-9.9s\t\t  -lnL0 = %.4f", model0->name, model0->ln);
    printf("\n   Alternative model = %-9.9s\t  -lnL1 = %.4f", model1->name, model1->ln);
    printf("\n   2(lnL1-lnL0) = %9.4f\t\t      df = %d ", delta, df);
    printf("\n   Using mixed chi-square distribution");
    if (prob == 1.0) {
        printf("\n   P-value = >%f", MAX_PROB);
    }
    else if (prob < 0.000001) {
        printf("\n   P-value = <%f", MIN_PROB);
    }
    else {
        printf("\n   P-value =  %f", prob);
    }
    return prob;
}

/******************* TestEqualBaseFrequencies ****************/
float TestEqualBaseFrequencies(ModelSt *model0, ModelSt *model1)
{
    float P;

    printf("\n Equal base frequencies");
    P = LRT(model0, model1);
    return P;
}

/*******************  TestTiequalsTv  **********************/
float TestTiequalsTv(ModelSt *model0, ModelSt *model1)
{
    float P;

    printf("\n Ti=Tv");
    P = LRT(model0, model1);
    return P;
}

/******************* TestEqualTiAndEqualTvRates *******************/
float TestEqualTiAndEqualTvRates (ModelSt *model0, ModelSt *model1)
{
    float P;

    printf("\n Unequal Tv and unequal Ti");
    P = LRT(model0, model1);
    return P;
}

/********************* TestEqualSiteRates **********************/
float TestEqualSiteRates (ModelSt *model0, ModelSt *model1)
{
    float P;

    printf("\n Equal rates among sites");
    if (mixchi) {
        P = LRTmix(model0, model1);
    }
    else {
        P = LRT(model0, model1);
    }
    return P;
}

/******************** TestInvariableSites **********************/
float TestInvariableSites (ModelSt *model0, ModelSt *model1)
{
    float P;

    printf("\n No Invariable sites");
    if (mixchi) {
        P = LRTmix(model0, model1);
    }
    else {
        P = LRT(model0, model1);
    }
    return P;
}

/**************  ChiSquare: probability of chi square value *************/
/*
ALGORITHM Compute probability of chi square value.
Adapted from:     Hill, I. D. and Pike, M. C.  Algorithm 299.Collected Algorithms for the CACM 1967 p. 243
Updated for rounding errors based on remark inACM TOMS June 1985, page 185. Found in Perlman.lib
*/
float ChiSquare (float x, int df)  /* x: obtained chi-square value,  df: degrees of freedom */
{
    float a, y, s;
    float e, c, z;
    int even;         /* true if df is an even number */

    if (x <= 0.0 || df < 1) {
        return (1.0);
    }
    y = 1;
    a = 0.5 * x;
    even = (2*(df/2)) == df;
    if (df > 1) {
        y = ex (-a);
    }
    s = (even ? y : (2.0 * Normalz (-sqrt(x))));
    if (df > 2) {
        x = 0.5 * (df - 1.0);
        z = (even ? 1.0 : 0.5);
        if (a > BIGX) {
            e = (even ? 0.0 : LOG_SQRT_PI);
            c = log (a);
            while (z <= x) {
                e = log (z) + e;
                s += ex (c*z-a-e);
                z += 1.0;
            }
            return (s);
        }
        else {
            e = (even ? 1.0 : (I_SQRT_PI / sqrt (a)));
            c = 0.0;
            while (z <= x) {
                e = e * (a / z);
                c = c + e;
                z += 1.0;
            }
            return (c * y + s);
        }
    }
    else {
        return (s);
    }
}

/************** Normalz: probability of normal z value *********************/
/*
ALGORITHM:    Adapted from a polynomial approximation in:
            Ibbetson D, Algorithm 209
            Collected Algorithms of the CACM 1963 p. 616
        Note:
            This routine has six digit accuracy, so it is only useful for absolute
            z values < 6.  For z values >= to 6.0, Normalz() returns 0.0.
*/
float Normalz (float z)        /*VAR returns cumulative probability from -oo to z VAR normal z value */
{
    float y, x, w;

    if (z == 0.0) {
        x = 0.0;
    }
    else {
        y = 0.5 * fabs (z);
        if (y >= (Z_MAX * 0.5)) {
            x = 1.0;
        }
        else if (y < 1.0) {
            w = y*y;
            x = ((((((((0.000124818987 * w
                -0.001075204047) * w +0.005198775019) * w
                -0.019198292004) * w +0.059054035642) * w
                -0.151968751364) * w +0.319152932694) * w
                -0.531923007300) * w +0.797884560593) * y * 2.0;
        }
        else {
            y -= 2.0;
            x = (((((((((((((-0.000045255659 * y
                +0.000152529290) * y -0.000019538132) * y
                -0.000676904986) * y +0.001390604284) * y
                -0.000794620820) * y -0.002034254874) * y
                +0.006549791214) * y -0.010557625006) * y
                +0.011630447319) * y -0.009279453341) * y
                +0.005353579108) * y -0.002141268741) * y
                +0.000535310849) * y +0.999936657524;
        }
    }
    return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}

/********************** RatioCalc ***************************/
static void RatioCalc()
{
    float score1, score2, ratio, prob;
    int df;

    printf("\nPlease, input the POSITIVE log likelihood score corresponding to \nthe null model> ");
    scanf("%f", &score1);
    while (score1 < 0) {
        printf("\nBad Input: the program doesn't accept negative likelihood scores");
        printf("\n\nPlease, input the POSITIVE log likelihood score corresponding to \nthe null model> ");
        scanf("%f", &score1);
    }
    printf("\nPlease, input the POSITIVE log likelihood score corresponding to \nthe alternative model> ");
    scanf("%f", &score2);
    while (score2 < 0) {
        printf("Bad Input: the program doesn't accept negative likelihood scores");
        printf("\nPlease, input the POSITIVE log likelihood score corresponding to \nthe alternative model> ");
        scanf("%f", &score2);
    }
    if (score1 < score2) {
        printf("\n\nIncorrect input: the positive likelihood of the null model cannot be smaller than the positive likelihood of the alternative model.");
        printf("\nYou should enter the correct positive log likelihood scores again.\n");
        exit(0);
    }
    printf("\nPlease, input the number of degrees of freedom> ");
    scanf("%d", &df);
    while (df < 1) {
        printf("\nThe number of degrees of freedom should be at least 1");
        printf("\nPlease, input the number of degrees of freedom> ");
        scanf("%d", &df);
    }
    ratio = 2*(score1-score2);
    prob = ChiSquare(ratio,df);
    printf("\n\n_________________________ Results of Ratio Calculator _______________________\n");
    printf("\nThe ratio is %f ", ratio);
    printf("\n\nThe probability of observing this ratio likelihood test statistic under a correct null model is %f\n", prob);
    if (prob < alpha) {
        printf ("\nThis is significant at the alpha level of %.4f", alpha);
    }
    else {
        printf ("\nThis is not significant at the alpha level of %.4f\n", alpha);
    }
    exit(0);
}

/*********************** CalculateAIC ***************************/
/* Calculates the AIC or AICc value for each likelihood score */
void CalculateAIC()
{
    float smallerAIC;
    int i, K, n;

    n = sampleSize;
    for (i = 0; i < NUM_MODELS; i++) {
        K = model[i].parameters;
        AIC[i] =  2 * (model[i].ln  +  K);
        if (useAICc == YES) {
            AIC[i] += 2*K*(K+1) / (double) (n-K-1);
        }
    }
    smallerAIC = AIC[0];
    for (i = 1; i < NUM_MODELS; i++) {
        if (AIC[i] < smallerAIC) {
            smallerAIC = AIC[i];
        }
    }
    minAIC = BIGNUMBER;
    for (i = NUM_MODELS - 1; i >= 0; i--) {
        if (AIC[i] == smallerAIC) {
            strcpy(modelAIC, model[i].name);
            minAIC = AIC[i];
        }
    }
}

/*********************** AkaikeWeights ****************************/
/* calculates deltaAIC and Akaike weights (w[i]) */
void AkaikeWeights ()
{
    int i, j, sorted, pass;
    float deltaAIC[NUM_MODELS];
    float ord[NUM_MODELS];
    float sumExp, temp;
    float cumWeight;

    sumExp = 0;
    for (i = 0; i < NUM_MODELS; i++) {
        deltaAIC[i] = AIC[i] - minAIC;
        sumExp += exp(-0.5 * deltaAIC[i]);
    }
    for (i = 0; i < NUM_MODELS; i++) {
        wAIC[i] = exp(-0.5 * deltaAIC[i]) / sumExp;
        ord[i] = AIC[i];
        orderedAIC[i] = i;
    }
    /* Sort by AIC score to print weights in order*/
    sorted = NO;
    pass = 1;
    while (sorted == NO) {
        sorted = YES;
        for (i = 0; i < (NUM_MODELS - pass); i++) {
             if (ord[i] > ord [i + 1]) {
                temp = ord[i + 1];
                ord[i + 1]= ord[i];
                ord[i] = temp;
                temp = orderedAIC[i + 1];
                orderedAIC[i + 1] = orderedAIC[i];
                orderedAIC[i] = temp;
                sorted = NO;
              }
        }
        pass++;
    }
    printf ("\n\n\n ** MODEL SELECTION UNCERTAINTY : Akaike Weights **");
    if (useAICc == NO) {
        printf ("\n\nModel\t\t-lnL\t\tK\t AIC\t\t delta\t\tWeight\t\tCumWeight");
    }
    else {
        printf ("\n\nModel\t\t-lnL\t\tK\t AICc\t\t delta\t\tWeight\t\tCumWeight");
    }
    printf ("\n-------------------------------------------------------------------------------------------------");
    cumWeight = 0;
    for (i = 0; i < NUM_MODELS; i++) {
        j = orderedAIC[i];
        cumWeight += wAIC[j];
        if (wAIC[j] > 0.0001) {
            printf("\n%-10s\t%10.4f\t%2d\t%10.4f\t%9.4f\t%8.4f\t%7.4f", model[j].name, model[j].ln, model[j].parameters, AIC[j], deltaAIC[j], wAIC[j], cumWeight);
        }
        else {
            printf("\n%-10s\t%10.4f\t%2d\t%10.4f\t%9.4f\t%4.2e\t%7.4f", model[j].name, model[j].ln, model[j].parameters, AIC[j], deltaAIC[j], wAIC[j], cumWeight);
        }
    }
    printf ("\n-------------------------------------------------------------------------------------------------");
    printf ("\n-lnL:\t\tnegative log likelihood");
    printf ("\n K:\t\tnumber of estimated (free) parameters");
    printf ("\n AIC:\t\tAkaike Information Criterion");
    printf ("\n delta:\t\tAkaike difference");
    printf ("\n weight:\tAkaike weight");
    printf ("\n cumWeight:\tcumulative Akaike weight");
}

// /************** ModelAveraging **********************/
// /*  Calculates the importance for different parameters
//     of the models (it is simply the sum of the Akaike
//     weights for those models that include such parameter)
//     and model averaged estimates
// 
//     This method is completely brute force (See Java version)
// 
//     Assumes TrN y TIM estimate only Rb, Re
//     K81 estimates no R parameter
//     TVM estimates only Ra, Rc, Rd
//     GTR y SIM estimate Ra, Rb, Rc, Rd, Re
// */
// 
// 
// void ModelAveraging()
// {
//     double minWeightToAverage;
//     double ifA, ifC, ifG, ifT, ititv, iRa, iRb, iRc, iRd, iRe, ipinvI, ialphaG, ipinvIG, ialphaIG;
//     double wfA, wfC, wfG, wfT, wtitv, wRa, wRb, wRc, wRd, wRe, wpinvI, walphaG, wpinvIG, walphaIG;
// 
//     /* which index (1-167) for scores */
//     int efA[] = {14,20,27,34,58,65,73,81,122,133,145,157};
//     int efC[] = {15,21,28,35,59,66,74,82,123,134,146,158};
//     int efG[] = {16,22,29,36,60,67,75,83,124,135,147,159};
//     int efT[] = {17,23,30,37,61,68,76,84,125,136,148,160};
//     int etitv[] = {42,45,49,53,62,69,77,85};
//     int eRa[] = {90,97,105,113,126,137,149,161};
//     int eRb[] = {91,98,106,114,127,138,150,162};
//     int eRc[] = {92,99,107,115,128,139,151,163};
//     int eRd[] = {93,100,108,116,129,140,152,164};
//     int eRe[] = {94,101,109,117,130,141,153,165};
//     int epinvI[] = {4,24,46,70,102,142};
//     int ealphaG[] = {7,31,50,78,110,154};
//     int epinvIG[] = {4,10,24,38,46,54,70,86,102,118,142,166};
//     int ealphaIG[] = {7,11,31,39,50,55,78,87,110,119,154,167};
// 
//     /* which index (1-23) for models containing the parameter  */
//     int mfA[] = {4,5,6,7,12,13,14,15,20,21,22,23};
//     int mfC[] = {4,5,6,7,12,13,14,15,20,21,22,23};
//     int mfG[] = {4,5,6,7,12,13,14,15,20,21,22,23};
//     int mfT[] = {4,5,6,7,12,13,14,15,20,21,22,23};
//     int mtitv[] = {8,9,10,11,12,13,14,15};
//     int mRa[] = {16,17,18,19,20,21,22,23};
//     int mRb[] = {16,17,18,19,20,21,22,23};
//     int mRc[] = {16,17,18,19,20,21,22,23};
//     int mRd[] = {16,17,18,19,20,21,22,23};
//     int mRe[] = {16,17,18,19,20,21,22,23};
//     int mpinvI[] = {1,5,9,13,17,21};
//     int malphaG[] = {2,6,10,14,18,22};
//     int mpinvIG[] = {1,3,5,7,9,11,13,15,17,19,21,23};
//     int malphaIG[] = {2,3,6,7,10,11,14,15,18,19,22,23};
// 
//     ifA = ifC = ifG = ifT = ititv = iRa = iRb = iRc = iRd = iRe = ipinvI = ialphaG = ipinvIG = ialphaIG = 0;
//     wfA = wfC = wfG = wfT = wtitv = wRa = wRb = wRc = wRd = wRe = wpinvI = walphaG = wpinvIG = walphaIG = 0;
// 
//     if (averagingConfidenceInterval < 1) {
//         minWeightToAverage = FindMinWeightToAverage ();
//     }
//     else {
//         minWeightToAverage = 0.0;
//         cumConfidenceWeight = 1.0;
//     }
// 
//     /* calculate importances and model-averaged estimates */
//     AverageEstimates ("fA",12, mfA, efA, &ifA, &wfA, minWeightToAverage);
//     AverageEstimates ("fC",12, mfC, efC, &ifC, &wfC, minWeightToAverage);
//     AverageEstimates ("fG",12, mfG, efG, &ifG, &wfG, minWeightToAverage);
//     AverageEstimates ("fT",12, mfT, efT, &ifT, &wfT, minWeightToAverage);
//     AverageEstimates ("titv", 8, mtitv, etitv, &ititv, &wtitv, minWeightToAverage);
//     AverageEstimates ("Ra",8, mRa, eRa, &iRa, &wRa, minWeightToAverage);
//     AverageEstimates ("Rb",8, mRb, eRb, &iRb, &wRb, minWeightToAverage);
//     AverageEstimates ("Rc",8, mRc, eRc, &iRc, &wRc, minWeightToAverage);
//     AverageEstimates ("Rd",8, mRd, eRd, &iRd, &wRd, minWeightToAverage);
//     AverageEstimates ("Re",8, mRe, eRe, &iRe, &wRe, minWeightToAverage);
//     AverageEstimates ("pinv(I)",6, mpinvI, epinvI, &ipinvI, &wpinvI, minWeightToAverage);
//     AverageEstimates ("alpha(G)",6, malphaG, ealphaG, &ialphaG, &walphaG, minWeightToAverage);
//     AverageEstimates ("pinv(IG)",12, mpinvIG, epinvIG, &ipinvIG, &wpinvIG, minWeightToAverage);
//     AverageEstimates ("alpha(IG)",12, malphaIG, ealphaIG, &ialphaIG, &walphaIG, minWeightToAverage);
// 
//     /* print results */
//     printf ("\n\n\n\n* MODEL AVERAGING AND PARAMETER IMPORTANCE (using Akaike Weights)");
//     if (averagingConfidenceInterval == 1) {
//         fprintf (stdout, "\n  Including all %d models", NUM_MODELS);
//     }
//     else {
//         fprintf (stdout, "\n    Including only the best %d models within the aproximate %4.2f (%6.4f)\n    confidence interval", lastModelConfidence+1, averagingConfidenceInterval, cumConfidenceWeight);
//         fprintf (stdout, "\n      minimum weight to average is %6.4f", minWeightToAverage);
//         fprintf (stdout, "\n      weights are reescaled by the interval cumulative weight (%6.4f)", cumConfidenceWeight);
//     }
// 
//     printf ("\n\n\t\t\t\t\tModel-averaged");    /*Att göra: fixa till output för printet nedan*/
//     printf ("\nParameter\t\tImportance\testimates");
//     printf ("\n----------------------------------------------------");
//     printf ("\nfA\t\t\t%6.4f\t\t%11s",ifA, CheckNA(wfA));
//     printf ("\nfC\t\t\t%6.4f\t\t%11s",ifC, CheckNA(wfC));
//     printf ("\nfG\t\t\t%6.4f\t\t%11s",ifG, CheckNA(wfG));
//     printf ("\nfT\t\t\t%6.4f\t\t%11s",ifT, CheckNA(wfT));
//     printf ("\nTiTv\t\t\t%6.4f\t\t%11s",ititv, CheckNA(wtitv));
//     printf ("\nrAC\t\t\t%6.4f\t\t%11s",iRa, CheckNA(wRa));
//     printf ("\nrAG\t\t\t%6.4f\t\t%11s",iRb, CheckNA(wRb));
//     printf ("\nrAT\t\t\t%6.4f\t\t%11s",iRc, CheckNA(wRc));
//     printf ("\nrCG\t\t\t%6.4f\t\t%11s",iRd, CheckNA(wRd));
//     printf ("\nrCT\t\t\t%6.4f\t\t%11s",iRe, CheckNA(wRe));
//     printf ("\npinv(I)\t\t\t%6.4f\t\t%11s",ipinvI, CheckNA(wpinvI));
//     printf ("\nalpha(G)\t\t%6.4f\t\t%11s",ialphaG, CheckNA(walphaG));
//     printf ("\npinv(I+IG)\t\t%6.4f\t\t%11s",ipinvIG, CheckNA(wpinvIG));
//     printf ("\nalpha(G+IG)\t\t%6.4f\t\t%11s",ialphaIG, CheckNA(walphaIG));
//     printf ("\n----------------------------------------------------");
// 
//     printf ("\nNote: values have been rounded.");
//     printf ("\n (I):\t\taveraged using only +I models");
//     printf ("\n (G):\t\taveraged using only +G models");
//     printf ("\n (I+IG):\taveraged using both +I and +I+G models");
//     printf ("\n (G+IG):\taveraged using both +G and +I+G models");
// 
// }


// /************** AverageEstimates **********************/
// /*
//     Calculates parameter importance and averaged estimates
// */
// void AverageEstimates (char *whichParameter, int numModels, int *modelIndex, int *estimateIndex, double *importance, double *averagedEstimate, double minWeightToAverage)
// {
//     int i;
// 
//     whichParameter = whichParameter; /* just to avoid warnings */
//     for (i=0; i < numModels; i++) {
//         if (wAIC[modelIndex[i]] < minWeightToAverage) {
//             continue;
//         }
//         *importance += wAIC[modelIndex[i]];
//         *averagedEstimate += wAIC[modelIndex[i]] * score[estimateIndex[i]];
//     }
//     /* rescale importance to the total weight of the models included in the confidence interval */
//     if (*importance  > 0) {
//         *averagedEstimate /= *importance;
//         *importance /= cumConfidenceWeight;
//     }
//     else {
//         *averagedEstimate = NA;
//     }
// }
// 
// /*********************** FindMinWeightToAverage ****************************/
// /*
//     Finds the minimum weight we want to average, that is the weight corresponding
//     to the last model in the confidence interval specified
// */
// double FindMinWeightToAverage ()
// {
//     int i;
//     double minWeight,cumWeight;
// 
//     cumWeight = 0;
//     for (i = 0; i < NUM_MODELS; i++) {
//         cumWeight += wAIC[orderedAIC[i]];
//         if (cumWeight > averagingConfidenceInterval) {
//             minWeight = wAIC[orderedAIC[i]];
//             lastModelConfidence = i;
//             cumConfidenceWeight = cumWeight;
//             break;
//         }
//     }
//     return minWeight;
// }

/*********************** AICfile ****************************/
/* reads likelihood scores from a file and calculates their AIC */
/* values, choosing the minimum */
void AICfile()
{
    float ln[100];
    int n[100];
    float AIC[100];
    float min_AIC;
    int i, j;

    i = 0;
    printf("\n AIC calculation from file \n");
    while (!feof(stdin)) {
        i++;
        scanf("%f %d", &ln[i], &n[i]);
        AIC[i] = 2*(ln[i] + n[i]);
    }
    i--;
    min_AIC = AIC[1];
    for (j = 1; j <= i; j++) {
        if (AIC[j] <= min_AIC) {
            min_AIC = AIC[j];
        }
    }
    printf("\nNumber\t\tLikelihood\t\tParameters\t\tAIC\n");
    for (j = 1; j <= i; j++) {
        printf("%2d\t%15.5f\t%5d\t%15.5f\n", j, ln[j], n[j], AIC[j]);
    }
    for (j = 1; j <= i; j++) {
        if (AIC[j] == min_AIC) {
            printf("\n A minimum AIC value (%f) corresponds to the score number %d (%f)",  min_AIC, j, ln[j]);
        }
    }
    exit(0);
}

/*********************** AICCalc ***************************/
/*Ask the user for likelihood scores and calculates their AIC */
/*values, choosing the minimum. */
//void AICCalc()
//{
//    double ln[300];
//    int n[300];
//    double AIC[300];
//    double min_AIC;
//    int i, number;
//
//    printf("\nEnter the number of likelihood scores you want to compare> ");
//    scanf("%d", &number);
//    for (i = 1; i <= number; i++) {
//        printf("\nEnter the positive likelihood score number %d> ", i);
//        scanf("%lf", &ln[i]);
//        printf("\nEnter the number of free parameters corresponding to the model \nrepresented by the score number %d> ", i);
//        scanf("%d", &n[i]);
//        if (ln[i] < 0 || n[i] < 0 || ln[i] == 0 || n[i] == 0) {
//            printf("\nThe program only admits positive likelihood scores or number of parameters.");
//            printf("\nEnter the positive likelihood score number %d", i);
//            scanf("%lf", &ln[i]);
//        }
//        AIC[i] = 2*(ln[i] + n[i]);
//    }
//    min_AIC = AIC[1];
//    for (i = 1; i <= number; i++) {
//        if (AIC[i] <= min_AIC) {
//            min_AIC = AIC[i];
//        }
//    }
//    printf("\n\n_________________________ Results of AIC Calculator _______________________\n");
//    printf("\nNumber\t\tLikelihood\t\tParameters\t\tAIC\n");
//    for (i = 1; i <= number; i++) {
//        printf("%2d\t%15.5f\t%5d\t%15.5f\n", i, ln[i], n[i], AIC[i]);
//    }
//    for (i = 1; i <= number; i++) {
//        if (AIC[i] == min_AIC) {
//            printf("\n A minimum AIC value (%f) corresponds to the score number %d (%f)\n",  min_AIC, i, ln[i]);
//        }
//    }
//    printf("\nDone.\n");
//    exit(0);
//}
int AICCalc()
{

    double ln[300];
    int n[300];
    double AIC[300];
    int min_AIC;
    int i, number;

    while (1) {
        fprintf(stderr,
            "Enter the number of likelihood scores you want to compare> ");
        if (scanf("%d", &number) != 1) {
            perror("Error in scanf() call");
            exit(EXIT_FAILURE);
        }
        else if (number > 0 && number < 300) {
            break;
        }
        fprintf(stderr, "Invalid number, must be 1 or larger (but not as large as 300)\n");
    }
    for (i = 0; i < number; i++) {
        putchar('\n');
        fprintf(stderr, "Enter the positive likelihood score number %d> ", i + 1);
        if (scanf("%lf", &ln[i]) != 1) {
            perror("Error in scanf()");
            exit(EXIT_FAILURE);
        }
        fprintf(stderr,
            "Enter the number of free parameters corresponding to the model "
            "\nrepresented by the score number %d> ",
            i + 1);
        if (scanf("%d", &n[i]) != 1) {
            perror("Error in scanf()");
            exit(EXIT_FAILURE);
        }
        if (ln[i] <= 0 || n[i] <= 0) {
            fputs("The program only admits positive likelihood scores or number "
                "of parameters\n",
                stderr);
            fprintf(stderr, "Enter the positive likelihood score number %d", i + 1);
        if (scanf("%lf", &ln[i]) != 1) {
            perror("Error in scanf()");
            exit(EXIT_FAILURE);
            }
        }
        AIC[i] = 2 * (ln[i] + n[i]);
    }
    min_AIC = 0;
    for (i = 1; i < number; i++) {
        if (AIC[i] < AIC[min_AIC]) {
            min_AIC = i;
        }
    }
    puts("\n");
    puts("_________________________ "
        "Results of AIC Calculator "
        "_______________________");
    puts("\nNumber\t\tLikelihood\t\tParameters\t\tAIC");
    for (i = 0; i < number; i++) {
        printf("%2d\t%15.5f\t%5d\t%15.5f\n", i + 1, ln[i], n[i], AIC[i]);
    }
    printf("\n A minimum AIC value (%f) corresponds to the score number %d "
        "(%f)\n",
        AIC[min_AIC], min_AIC + 1, ln[min_AIC]);
    puts("\nDone.");

    return EXIT_SUCCESS;
}

/********************* PrintPaupBlock ************************/
/* Prints a block of paup commands for appending to the data file */
static void PrintPaupBlock (int ishLRT)
{
    fT = 1 - (fA + fC + fG);
    printf("\n\n\n--\n\nPAUP* Commands Block:");
    printf(" If you want to implement the previous estimates as likelihod settings in PAUP*,");
    printf(" attach the next block of commands after the data in your PAUP file:\n");
    if (ishLRT == YES) {
        printf("\n[!\nLikelihood settings from best-fit model (%s) selected by hLRT in %s %s\n]", modelhLRT, PROGRAM_NAME, VERSION_NUMBER);
    }
    else if (useAICc == NO) {
        printf("\n[!\nLikelihood settings from best-fit model (%s) selected by AIC in %s %s\n]", modelAIC, PROGRAM_NAME, VERSION_NUMBER);
    }
    else {
        printf("\n[!\nLikelihood settings from best-fit model (%s) selected by AICc in %s %s\n]", modelAIC, PROGRAM_NAME, VERSION_NUMBER);
    }
    printf("\nBEGIN PAUP;");
    printf("\n\tLset");
    printf("  Base=");
    if (fA == fC && fA ==fG && fA == fT) {
        printf("equal");
    }
    else {
        printf("(%.4f %.4f %.4f)",fA,fC,fG);
    }
    /* Substitution rates */
    /*if (Ra == Rb && Ra == Rc && Ra == Rd && Ra == Re && Ra == Rf && TiTv == 0) {*/
    if (rAC == rAG && rAC == rAT && rAC == rCG && rAC == rCT && rAC == rGT && TiTv == 0) {
        printf("  Nst=1");
    }
    else if (TiTv != 0) {
        printf("  Nst=2  TRatio=%.4f", TiTv);
    }
    else {
        /*printf("  Nst=6  Rmat=(%.9f %.9f %.9f %.9f %.9f)", Ra, Rb, Rc, Rd, Re);*/
        printf("  Nst=6  Rmat=(%.9f %.9f %.9f %.9f %.9f)", rAC, rAG, rAT, rCG, rCT);
    }
    /* Rate variation */
    printf("  Rates=");
    if (shape == 0 || shape > 999) {
        printf("equal");
    }
    else {
        printf("gamma  Shape=%.4f", shape);
    }
    /* Invariable sites */
    printf("  Pinvar=");
    if (pinv == 0) {
        printf("0");
    }
    else {
        printf("%.4f", pinv);
    }
    printf(";\nEND;");
    printf("\n\n--");
}

/********************* PrintMbBlock ************************/
/* Prints a block of MrBayes commands for appending to the data file */
static void PrintMbBlock (int ishLRT)
{
    fT = 1 - (fA + fC + fG);
    printf("\n\n\nMrBayes Commands Block:");
    printf(" If you want to implement a \"best\" model in MrBayes,");
    printf(" attach the next block of commands after the data in your NEXUS file:\n");
    printf("(NOTE: In a Bayesian analysis, the Markov chain is integrating over the");
    printf(" uncertainty in parameter values. Thus, you usually do NOT want to use");
    printf(" the parameter values estimated by the commands in MrModeltest or Modeltest.");
    printf(" You rather want to specify the general \"form\" of the model (such as nst=1 etc.)\n");
    if (ishLRT == YES) {
        printf("\n[!\nMrBayes settings for the best-fit model (%s) selected by hLRT in %s %s\n]", modelhLRT, PROGRAM_NAME, VERSION_NUMBER);
    }
    else if (useAICc == NO) {
        printf("\n[!\nMrBayes settings for the best-fit model (%s) selected by AIC in %s %s\n]", modelAIC, PROGRAM_NAME, VERSION_NUMBER);
    }
    else {
        printf("\n[!\nMrBayes settings for the best-fit model (%s) selected by AICc in %s %s\n]", modelAIC, PROGRAM_NAME, VERSION_NUMBER);
    }
    printf("\nBEGIN MRBAYES;\n");
    printf("\n\tLset");
    /* Substitution rates */
    /*if (Ra == Rb && Ra == Rc && Ra == Rd && Ra == Re && Ra == Rf && TiTv == 0) {*/
    if (rAC == rAG && rAC == rAT && rAC == rCG && rAC == rCT && rAC == rGT && TiTv == 0) {
        printf("  nst=1");
    }
    else if (TiTv != 0) { /*** Att göra: Monitor this. Might be "yes" also for nst=6! **/
        printf("  nst=2");
    }
    else {
        printf("  nst=6");
    }
    /* Rate variation */
    printf("  rates=");
    if (pinv == 0) {
        if (shape == 0 || shape > 999) {
            printf("equal");
        }
        else {
            printf("gamma");
        }
    }
    else if (shape == 0 || shape > 999) {
        printf("propinv");
    }
    else {
        printf("invgamma");
    }
    printf(";\n");
    /* Base frequencies */
    if (fA == fC && fA == fG && fA == fT) {
        printf("\tPrset statefreqpr=fixed(equal);");
    }
    else {
        printf("\tPrset statefreqpr=dirichlet(1,1,1,1);");
    }
    printf("\nEND;");
    printf("\n\n--");
}

/********************* HLRTAttention ************************/
/* [Not completed] Warn if the different hLRT hierarchies give different models  */
static void HLRTAttention(char *first, char *second, char *third, char *fourth)
{
        printf("\n\n\n --");
        printf("\n ATTENTION: The choice based on hLRT can be sensitive for the specific");
        printf("\n            hierarchy used. If selected models differ, User need to");
        printf("\n            make the choice!");
        printf("\n\n          Model selected by hLRT (default): %s", first);
        printf("\n            Model selected by hLRT2:          %s", second);
        printf("\n            Model selected by hLRT3:          %s", third);
        printf("\n            Model selected by hLRT4:          %s", fourth);
        printf("\n --\n");
}

/********************* Output ************************/
/* Prints the results of MrModeltest  */
static void Output(char *selection, float value)
{
    int i, numK;

    fT = 1 - (fA + fG + fC);
    for (i = 0; i < NUM_MODELS; i++) {
        if (!strcmp (selection, model[i].name)) {
            theln = model[i].ln;
            numK = model[i].parameters;
        }
    }
    printf("\n\n Model selected: %s", selection);
    printf("\n   -lnL = \t%7.4f", theln);
    printf("\n    K = \t%d", numK);
    if (value > 0) {
        if (useAICc == YES) {
            printf("\n    AICc = \t%7.4f\n", value);
        }
        else {
            printf("\n    AIC = \t%7.4f\n", value);
        }
    }
    printf("\n   Base frequencies: ");
    if (fA == fC && fA == fG && fA == fT) {
        printf("\n     Equal frequencies");
    }
    else {
        printf("\n     freqA = \t%7.4f", fA);
        printf("\n     freqC = \t%7.4f", fC);
        printf("\n     freqG = \t%7.4f", fG);
        printf("\n     freqT = \t%7.4f", fT);
    }
    printf("\n   Substitution model: ");
    /*if (Ra == Rb && Ra == Rc && Ra == Rd && Ra == Re && Ra == Rf && TiTv == 0) {*/
    if (rAC == rAG && rAC == rAT && rAC == rCG && rAC == rCT && rAC == rGT && TiTv == 0) {
        printf("\n     All rates equal");
    }
    else if (TiTv != 0) {
        printf("\n    Ti/tv ratio =\t%7.4f", TiTv);
    }
    else {
        printf("\n     Rate matrix");
        /*printf("\n     R(a) [A-C] = \t%7.4f", Ra);*/
        /*printf("\n     R(b) [A-G] = \t%7.4f", Rb);*/
        /*printf("\n     R(c) [A-T] = \t%7.4f", Rc);*/
        /*printf("\n     R(d) [C-G] = \t%7.4f", Rd);*/
        /*printf("\n     R(e) [C-T] = \t%7.4f", Re);*/
        /*printf("\n     R(f) [G-T] = \t%7.4f", 1.0)*/;
        printf("\n     rAC = \t%7.4f", rAC);
        printf("\n     rAG = \t%7.4f", rAG);
        printf("\n     rAT = \t%7.4f", rAT);
        printf("\n     rCG = \t%7.4f", rCG);
        printf("\n     rCT = \t%7.4f", rCT);
        printf("\n     rGT = \t%7.4f", rGT);
    }
    printf("\n   Among-site rate variation");
    if (pinv == 0) {
        printf("\n     Proportion of invariable sites = 0");
    }
    else {
        printf("\n     Proportion of invariable sites (I) = %.4f", pinv);
        printf("\n     Variable sites (G)");
    }
    if (shape == 0) {
        printf("\n     Equal rates for all sites");
    }
    else if (shape > 999) { /* shape is infinity */
        printf("\n     Equal rates for all sites (shape parameter = infinity)");
    }
    else {
        printf("\n     Gamma distribution shape parameter = %.4f", shape);
    }
}

/********************* SetModel ************************/
/* Sets the parameter estimates for the selected model */
static void SetModel(char *selection)
{
    /* Default parameter estimates for the selected model (JC)*/
    fA = fC = fG = fT = 0.25;
    TiTv = 0;
    /*Ra = Rb = Rc = Rd = Re = Rf = 1.0;*/
    rAC = rAG = rAT = rCG = rCT = rGT = 1.0;
    shape = 0.0;
    pinv = 0.0;
    if (!strcmp (selection, "JCI")) {
        pinv = score[4];
    }
    else if (!strcmp (selection, "JC+G")) {
        shape = score[7];
    }
    else if (!strcmp (selection, "JC+I+G")) {
        pinv = score[10];
        shape = score[11];
    }
    else if (!strcmp (selection, "F81")) {
        fA = score[14];
        fC = score[15];
        fG = score[16];
        fT = score[17];
    }
    else if (!strcmp (selection, "F81+I")) {
        fA = score[20];
        fC = score[21];
        fG = score[22];
        fT = score[23];
        pinv = score[24];
    }
    else if (!strcmp (selection, "F81+G")) {
        fA = score[27];
        fC = score[28];
        fG = score[29];
        fT = score[30];
        shape = score[31];
    }
    else if (!strcmp (selection, "F81+I+G")) {
        fA = score[34];
        fC = score[35];
        fG = score[36];
        fT = score[37];
        pinv = score[38];
        shape = score[39];
    }
    else if (!strcmp (selection, "K80")) {
        TiTv = score[42];
    }
    else if (!strcmp (selection, "K80+I")) {
        TiTv = score[45];
        pinv = score[46];
    }
    else if (!strcmp (selection, "K80+G")) {
        TiTv = score[49];
        shape = score[50];
    }
    else if (!strcmp (selection, "K80+I+G")) {
        TiTv = score[53];
        pinv = score[54];
        shape = score[55];
    }
    else if (!strcmp (selection, "HKY")) {
        fA = score[58];
        fC = score[59];
        fG = score[60];
        fT = score[61];
        TiTv = score[62];
    }
    else if (!strcmp (selection, "HKY+I")) {
        fA = score[65];
        fC = score[66];
        fG = score[67];
        fT = score[68];
        TiTv = score[69];
        pinv = score[70];
    }
    else if (!strcmp (selection, "HKY+G")) {
        fA = score[73];
        fC = score[74];
        fG = score[75];
        fT = score[76];
        TiTv = score[77];
        shape = score[78];
        }
    else if (!strcmp (selection, "HKY+I+G")) {
        fA = score[81];
        fC = score[82];
        fG = score[83];
        fT = score[84];
        TiTv = score[85];
        pinv = score[86];
        shape = score[87];
        }
    else if (!strcmp (selection, "SYM")) {
        rAC = score[90]; /*Ra = score[90];*/
        rAG = score[91]; /*Rb = score[91];*/
        rAT = score[92]; /*Rc = score[92];*/
        rCG = score[93]; /*Rd = score[93];*/
        rCT = score[94]; /*Re = score[94];*/
        rGT = score[95];
    }
    else if (!strcmp (selection, "SYM+I")) {
        /*Ra = score[97];*/
        /*Rb = score[98];*/
        /*Rc = score[99];*/
        /*Rd = score[100];*/
        /*Re = score[101];*/
        rAC = score[98];
        rAG = score[99];
        rAT = score[100];
        rCG = score[101];
        rCT = score[102];
        rGT = score[103];
        pinv = score[104]; /*102*/
    }
    else if (!strcmp (selection, "SYM+G")) {
        /*Ra = score[105];*/
        /*Rb = score[106];*/
        /*Rc = score[107];*/
        /*Rd = score[108];*/
        /*Re = score[109];*/
        rAC = score[107];
        rAG = score[108];
        rAT = score[109];
        rCG = score[110];
        rCT = score[111];
        rGT = score[112];
        shape = score[113]; /*110*/
        }
    else if (!strcmp (selection, "SYM+I+G")) {
        /*Ra = score[113];*/
        /*Rb = score[114];*/
        /*Rc = score[115];*/
        /*Rd = score[116];*/
        /*Re = score[117];*/
        rAC = score[116];
        rAG = score[117];
        rAT = score[118];
        rCG = score[119];
        rCT = score[120];
        rGT = score[121];
        pinv = score[122]; /*118*/
        shape = score[123]; /*119*/
    }
    else if (!strcmp (selection, "GTR")) {
        fA = score[126]; /*122*/
        fC = score[127]; /*123*/
        fG = score[128]; /*124*/
        fT = score[129]; /*125*/
        /*Ra = score[126];*/
        /*Rb = score[127];*/
        /*Rc = score[128];*/
        /*Rd = score[129];*/
        /*Re = score[130];*/
        rAC = score[130];
        rAG = score[131];
        rAT = score[132];
        rCG = score[133];
        rCT = score[134];
        rGT = score[135];
    }
    else if (!strcmp (selection, "GTR+I")) {
        fA = score[138]; /*133*/
        fC = score[139]; /*134*/
        fG = score[140]; /*135*/
        fT = score[141]; /*136*/
        /*Ra = score[137];*/
        /*Rb = score[138];*/
        /*Rc = score[139];*/
        /*Rd = score[140];*/
        /*Re = score[141];*/
        rAC = score[142];
        rAG = score[143];
        rAT = score[144];
        rCG = score[145];
        rCT = score[146];
        rGT = score[147];
        pinv = score[148]; /*142*/
    }
    else if (!strcmp (selection, "GTR+G")) {
        fA = score[151]; /*145*/
        fC = score[152]; /*146*/
        fG = score[153]; /*147*/
        fT = score[154]; /*148*/
        /*Ra = score[149];*/
        /*Rb = score[150];*/
        /*Rc = score[151];*/
        /*Rd = score[152];*/
        /*Re = score[153];*/
        rAC = score[155];
        rAG = score[156];
        rAT = score[157];
        rCG = score[158];
        rCT = score[159];
        rGT = score[160];
        shape = score[161]; /*154*/
    }
    else if (!strcmp (selection, "GTR+I+G")) {
        fA = score[164]; /*157*/
        fC = score[165]; /*158*/
        fG = score[166]; /*159*/
        fT = score[167]; /*160*/
        rAC = score[168];
        rAG = score[169];
        rAT = score[170];
        rCG = score[171];
        rCT = score[172];
        rGT = score[173];
        /*Ra = score[161];*/
        /*Rb = score[162];*/
        /*Rc = score[163];*/
        /*Rd = score[164];*/
        /*Re = score[165];*/
        pinv = score[174]; /*166*/
        shape = score[175]; /*167*/
    }
}

/********************* hLRT  ************************/
/* Performs the hypothesis testing using the hLRT1 hierarchy in Posada & Crandall 2001 and print the results*/
void hLRT()
{
    if (TestEqualBaseFrequencies (JC, F81) < alpha) { /* 1,2 */
        if (TestTiequalsTv (F81, HKY) < alpha) { /* 3,4 */
            if (TestEqualTiAndEqualTvRates (HKY,GTR) < alpha) { /* 5,6 */
                if (TestEqualSiteRates (GTR, GTRG) < alpha) { /* 7, 8 */
                    if (TestInvariableSites (GTRG, GTRIG) < alpha) { /*  9,10  */
                        strcpy(modelhLRT,"GTR+I+G"); /* 12 */
                    }
                    else {
                        strcpy(modelhLRT,"GTR+G"); /* 11 */
                    }
                }
                else {
                    if (TestInvariableSites (GTR, GTRI) < alpha) { /* 13 14,  */
                        strcpy(modelhLRT,"GTR+I"); /* 16 */
                    }
                    else {
                        strcpy(modelhLRT,"GTR"); /* 15 */
                    }
                }
            }
            else {
                if (TestEqualSiteRates (HKY, HKYG) < alpha) { /* 17 , 18 */
                    if (TestInvariableSites (HKYG, HKYIG) < alpha) { /* 19 , 20 */
                        strcpy(modelhLRT,"HKY+I+G"); /* 22 */
                    }
                    else {
                        strcpy(modelhLRT,"HKY+G"); /* 21 */
                    }
                }
                else {
                    if (TestInvariableSites (HKY, HKYI) < alpha) { /*  23, 24 */
                        strcpy(modelhLRT,"HKY+I"); /* 26 */
                    }
                    else {
                        strcpy(modelhLRT,"HKY"); /* 25 */
                    }
                }
            }
        }
        else {
            if (TestEqualSiteRates (F81, F81G) < alpha) { /* 27 , 28 */
                if (TestInvariableSites (F81G, F81IG) < alpha) { /* 29 , 30 */
                        strcpy(modelhLRT,"F81+I+G"); /* 32 */
                }
                else {
                    strcpy(modelhLRT,"F81+G"); /* 31 */
                }
            }
            else {
                if (TestInvariableSites (F81, F81I) < alpha) { /* 33 , 34 */
                    strcpy(modelhLRT,"F81+I"); /* 36 */
                }
                else {
                    strcpy(modelhLRT,"F81"); /* 35 */
                }
            }
        }
    }
    else {
        if (TestTiequalsTv (JC, K80) < alpha) { /* 37 , 38 */
            if (TestEqualTiAndEqualTvRates (K80, SYM) < alpha) { /* 39 , 40 */
                if (TestEqualSiteRates (SYM, SYMG) < alpha) { /* 41 , 42 */
                    if (TestInvariableSites (SYMG, SYMIG) < alpha) { /*  43,  44*/
                        strcpy(modelhLRT,"SYM+I+G"); /* 46 */
                    }
                    else {
                        strcpy(modelhLRT,"SYM+G"); /* 45 */
                    }
                }
                else {
                    if (TestInvariableSites (SYM, SYMI) < alpha) { /* 47 , 48 */
                        strcpy(modelhLRT,"SYM+I"); /* 50 */
                    }
                    else {
                        strcpy(modelhLRT,"SYM"); /* 49 */
                    }
                }
            }
            else {
                if (TestEqualSiteRates (K80, K80G) < alpha) { /*  51, 52 */
                    if (TestInvariableSites (K80G, K80IG) < alpha) { /* 53 , 54 */
                        strcpy(modelhLRT,"K80+I+G"); /* 56 */
                    }
                    else {
                        strcpy(modelhLRT,"K80+G"); /* 55 */
                    }
                }
                else {
                    if (TestInvariableSites (K80, K80I) < alpha) { /*  57, 58 */
                        strcpy(modelhLRT,"K80+I"); /* 60 */
                    }
                    else {
                        strcpy(modelhLRT,"K80"); /* 59 */
                    }
                }
            }
        }
        else {
            if (TestEqualSiteRates (JC, JCG) < alpha) { /* 61 ,  62*/
                if (TestInvariableSites (JCG, JCIG) < alpha) { /* 63 , 64 */
                    strcpy(modelhLRT,"JC+I+G"); /* 66 */
                }
                else {
                    strcpy(modelhLRT,"JC+G"); /* 65 */
                }
            }
            else {
                if (TestInvariableSites (JC, JCI) < alpha) {/* 67, 68 */
                    strcpy(modelhLRT,"JC+I"); /* 70 */
                }
                else {
                    strcpy(modelhLRT,"JC"); /* 69 */
                }
            }
        }
    }
} /* end of method */

/********************* hLRT2 ************************/
/* Performs the hypothesis testing using the hLRT2 hierarchy in Posada & Crandall 2001. */
void hLRT2()
{
    if (TestEqualBaseFrequencies (SYMIG, GTRIG) < alpha) { /*A*/
        if (TestEqualTiAndEqualTvRates (HKYIG, GTRIG) < alpha) { /*B*/
            if (TestEqualSiteRates (GTRI, GTRIG) < alpha) { /*C*/
                if (TestInvariableSites (GTRG, GTRIG) < alpha) { /*D*/
                    strcpy(modelhLRT2,"GTR+I+G");
                }
                else {
                    strcpy(modelhLRT2,"GTR+G");
                }
            }
            else {
                if (TestInvariableSites (GTR, GTRI) < alpha) { /*E*/
                    strcpy(modelhLRT2,"GTR+I");
                }
                else {
                    strcpy(modelhLRT2,"GTR");
                }
            }
        }
        else {
            if (TestTiequalsTv (F81IG, HKYIG) < alpha) { /*F*/
                if (TestEqualSiteRates (HKYI, HKYIG) < alpha) { /*G*/
                    if (TestInvariableSites (HKYG, HKYIG) < alpha) { /*H*/
                        strcpy(modelhLRT2,"HKY+I+G");
                    }
                    else {
                        strcpy(modelhLRT2,"HKY+G");
                    }
                }
                else {
                    if (TestInvariableSites (HKY, HKYI) < alpha) { /*I*/
                        strcpy(modelhLRT2,"HKY+I");
                    }
                    else {
                        strcpy(modelhLRT2,"HKY");
                    }
                }
            }
            else {
                if (TestEqualSiteRates (F81I, F81IG) < alpha) { /*J*/
                    if (TestInvariableSites (F81G, F81IG) < alpha) { /*K*/
                        strcpy(modelhLRT2,"F81+I+G");
                    }
                    else {
                        strcpy(modelhLRT2,"F81+G");
                    }
                }
                else {
                    if (TestInvariableSites (F81, F81I) < alpha) { /*L*/
                        strcpy(modelhLRT2,"F81+I");
                    }
                    else {
                        strcpy(modelhLRT2,"F81");
                    }
                }
            }
        }
    }
    else {
        if (TestEqualTiAndEqualTvRates (K80IG, SYMIG) < alpha) { /*M*/
            if (TestEqualSiteRates (SYMI, SYMIG) < alpha) { /*N*/
                if (TestInvariableSites (SYMG, SYMIG) < alpha) { /*O*/
                    strcpy(modelhLRT2,"SYM+I+G");
                }
                else {
                    strcpy(modelhLRT2,"SYM+G");
                }
            }
            else {
                if (TestInvariableSites (SYM, SYMI) < alpha) { /*P*/
                    strcpy(modelhLRT2,"SYM+I");
                }
                else {
                    strcpy(modelhLRT2,"SYM");
                }
            }
        }
        else {
            if (TestTiequalsTv (JCIG, K80IG) < alpha) { /*Q*/
                if (TestEqualSiteRates (K80I, K80IG) < alpha) { /*R*/
                    if (TestInvariableSites (K80G, K80IG) < alpha) { /*S*/
                        strcpy(modelhLRT2,"K80+I+G");
                    }
                    else {
                        strcpy(modelhLRT2,"K80+G");
                    }
                }
                else {
                    if (TestInvariableSites (K80, K80I) < alpha) { /*T*/
                        strcpy(modelhLRT2,"K80+I");
                    }
                    else {
                        strcpy(modelhLRT2,"K80");
                    }
                }
            }
            else {
                if (TestEqualSiteRates (JCI, JCIG) < alpha) { /*U*/
                    if (TestInvariableSites (JCG, JCIG) < alpha) { /*V*/
                        strcpy(modelhLRT2,"JC+I+G");
                    }
                    else {
                        strcpy(modelhLRT2,"JC+G");
                    }
                }
                else {
                    if (TestInvariableSites (JC, JCI) < alpha) { /*X*/
                        strcpy(modelhLRT2,"JC+I");
                    }
                    else {
                        strcpy(modelhLRT2,"JC");
                    }
                }
            }
        }
    }
}

/********************* hLRT3 ************************/
/* Performs the hypothesis testing using the hLRT3 hierarchy in Posada & Crandall 2001. */
void hLRT3()
{
    if (TestEqualSiteRates (JC, JCG) < alpha) { /*A*/
        if (TestInvariableSites (JCG, JCIG) < alpha) { /*B*/
            if (TestTiequalsTv (JCIG, K80IG) < alpha) { /*C*/
                if (TestEqualTiAndEqualTvRates (K80IG, SYMIG) < alpha) { /*D*/
                    if (TestEqualBaseFrequencies (SYMIG, GTRIG) < alpha) { /*E*/
                        strcpy(modelhLRT3,"GTR+I+G");
                    }
                    else {
                        strcpy(modelhLRT3,"SYM+I+G");
                    }
                }
                else {
                    if (TestEqualBaseFrequencies (K80IG, HKYIG) < alpha) { /*F*/
                        strcpy(modelhLRT3,"HKY+I+G");
                    }
                    else {
                        strcpy(modelhLRT3,"K80+I+G");
                    }
                }
            }
            else {
                if (TestEqualBaseFrequencies (JCIG, F81IG) < alpha) { /*G*/
                    strcpy(modelhLRT3,"F81+I+G");
                }
                else {
                    strcpy(modelhLRT3,"JC+I+G");
                }
            }
        }
        else {
            if (TestTiequalsTv (JCG, K80G) < alpha) { /*H*/
                if (TestEqualTiAndEqualTvRates (K80G, SYMG) < alpha) { /*I*/
                    if (TestEqualBaseFrequencies (SYMG, GTRG) < alpha) { /*J*/
                        strcpy(modelhLRT3,"GTR+G");
                    }
                    else {
                        strcpy(modelhLRT3,"SYM+G");
                    }
                }
                else {
                    if (TestEqualBaseFrequencies (K80G, HKYG) < alpha) { /*K*/
                        strcpy(modelhLRT3,"HKY+G");
                    }
                    else {
                        strcpy(modelhLRT3,"K80+G");
                    }
                }
            }
            else {
                if (TestInvariableSites (JCG, F81G) < alpha) { /*L*/
                    strcpy(modelhLRT3,"F81+G");
                }
                else {
                    strcpy(modelhLRT3,"JC+G");
                }
            }
        }
    }
    else {
        if (TestInvariableSites (JC, JCI) < alpha) { /*M*/
            if (TestTiequalsTv (JCI, K80I) < alpha) { /*N*/
                if (TestEqualTiAndEqualTvRates (K80I, SYMI) < alpha) { /*O*/
                    if (TestEqualBaseFrequencies (SYMI, GTRI) < alpha) { /*P*/
                        strcpy(modelhLRT3,"GTR+I");
                    }
                    else {
                        strcpy(modelhLRT3,"SYM+I");
                    }
                }
                else {
                    if (TestEqualBaseFrequencies (K80I, HKYI) < alpha) { /*Q*/
                        strcpy(modelhLRT3,"HKY+I");
                    }
                    else {
                        strcpy(modelhLRT3,"K80+I");
                    }
                }
            }
            else {
                if (TestEqualBaseFrequencies (JCI, F81I) < alpha) { /*R*/
                    strcpy(modelhLRT3,"F81+I");
                }
                else {
                    strcpy(modelhLRT3,"JC+I");
                }
            }
        }
        else {
            if (TestTiequalsTv (JC, K80) < alpha) { /*S*/
                if (TestEqualTiAndEqualTvRates (K80, SYM) < alpha) { /*T*/
                    if (TestEqualBaseFrequencies (SYM, GTR) < alpha) { /*U*/
                        strcpy(modelhLRT3,"GTR");
                    }
                    else {
                        strcpy(modelhLRT3,"SYM");
                    }
                }
                else {
                    if (TestEqualBaseFrequencies (K80, HKY) < alpha) { /*V*/
                        strcpy(modelhLRT3,"HKY");
                    }
                    else {
                        strcpy(modelhLRT3,"K80");
                    }
                }
            }
            else {
                if (TestEqualBaseFrequencies (JC, F81) < alpha) { /*X*/
                    strcpy(modelhLRT3,"F81");
                }
                else {
                    strcpy(modelhLRT3,"JC");
                }
            }
        }
    }
}

/********************* hLRT4 ************************/
/* Performs the hypothesis testing using the hLRT4 hierarchy in Posada & Crandall 2001. */
void hLRT4()
{
    if (TestEqualSiteRates (GTRI, GTRIG) < alpha) { /*A*/
        if (TestInvariableSites (GTRG, GTRIG) < alpha) { /*B*/
            if (TestEqualTiAndEqualTvRates (HKYIG, GTRIG) < alpha) { /*C*/
                if (TestEqualBaseFrequencies (SYMIG, GTRIG) < alpha) { /*D*/
                    strcpy(modelhLRT4,"GTR+I+G");
                }
                else {
                    strcpy(modelhLRT4,"SYM+I+G");
                }
            }
            else {
                if (TestTiequalsTv (F81IG, HKYIG) < alpha) { /*E*/
                    if (TestEqualBaseFrequencies (K80IG, HKYIG) < alpha) { /*F*/
                        strcpy(modelhLRT4,"HKY+I+G");
                    }
                    else {
                        strcpy(modelhLRT4,"K80+I+G");
                    }
                }
                else {
                    if (TestEqualBaseFrequencies (JCIG, F81IG) < alpha) { /*G*/
                        strcpy(modelhLRT4,"F81+I+G");
                    }
                    else {
                        strcpy(modelhLRT4,"JC+I+G");
                    }
                }
            }
        }
        else {
            if (TestEqualTiAndEqualTvRates (HKYG, GTRG) < alpha) { /*H*/
                if (TestEqualBaseFrequencies (SYMG, GTRG) < alpha) { /*I*/
                    strcpy(modelhLRT4,"GTR+G");
                }
                else {
                    strcpy(modelhLRT4,"SYM+G");
                }
            }
            else {
                if (TestTiequalsTv (F81G, HKYG) < alpha) { /*J*/
                    if (TestEqualBaseFrequencies (K80G, HKYG) < alpha) { /*K*/
                        strcpy(modelhLRT4,"HKY+G");
                    }
                    else {
                        strcpy(modelhLRT4,"K80+G");
                    }
                }
                else {
                    if (TestEqualBaseFrequencies (JCG, F81G) < alpha) { /*L*/
                        strcpy(modelhLRT4,"F81+G");
                    }
                    else {
                        strcpy(modelhLRT4,"JC+G");
                    }
                }
            }
        }
    }
    else {
        if (TestInvariableSites (GTR, GTRI) < alpha) { /*M*/
            if (TestEqualTiAndEqualTvRates (HKYI, GTRI) < alpha) { /*N*/
                if (TestEqualBaseFrequencies (SYMI, GTRI) < alpha) { /*O*/
                    strcpy(modelhLRT4,"GTR+I");
                }
                else {
                    strcpy(modelhLRT4,"SYM+I");
                }
            }
            else {
                if (TestTiequalsTv (F81I, HKYI) < alpha) { /*P*/
                    if (TestEqualBaseFrequencies (K80I, HKYI) < alpha) { /*Q*/
                        strcpy(modelhLRT4,"HKY+I");
                    }
                    else {
                        strcpy(modelhLRT4,"K80+I");
                    }
                }
                else {
                    if (TestEqualBaseFrequencies (JCI, F81I) < alpha) { /*R*/
                        strcpy(modelhLRT4,"F81+I");
                    }
                    else {
                        strcpy(modelhLRT4,"JC+I");
                    }
                }
            }
        }
        else {
            if (TestEqualTiAndEqualTvRates (HKY, GTR) < alpha) { /*S*/
                if (TestEqualBaseFrequencies (SYM, GTR) < alpha) { /*T*/
                    strcpy(modelhLRT4,"GTR");
                }
                else {
                    strcpy(modelhLRT4,"SYM");
                }
            }
            else {
                if (TestTiequalsTv (F81, HKY) < alpha) { /*U*/
                    if (TestEqualBaseFrequencies (K80, HKY) < alpha) { /*V*/
                        strcpy(modelhLRT4,"HKY");
                    }
                    else {
                        strcpy(modelhLRT4,"K80");
                    }
                }
                else {
                    if (TestEqualBaseFrequencies (JC, F81) < alpha) { /*X*/
                        strcpy(modelhLRT4,"F81");
                    }
                    else {
                        strcpy(modelhLRT4,"JC");
                    }
                }
            }
        }
    }
}


/********************* PrintTitle **********************/
static void PrintTitle (FILE *fp)
{
    fprintf(fp, "\nOutput from %s version %s ", PROGRAM_NAME, VERSION_NUMBER);
    fprintf(fp, "\n\n%s is written by Johan Nylander and is a modified version of", PROGRAM_NAME);
    fprintf(fp, " Modeltest version 3.6 (Copyright David Posada, Universidad de Vigo).");
    fprintf(fp, "\n\nReference:");
    fprintf(fp, "\n\"Nylander, J.A.A. 2004. MrModeltest %s. Program", VERSION_NUMBER);
    fprintf(fp, " distributed by the author. Evolutionary Biology Centre, Uppsala University.\"");
    fprintf(fp, "\n\nContact:");
    fprintf(fp, "\njohan.nylander@ebc.uu.se.");
    fprintf(fp, "\n\nCredits: David Posada is thanked for supplying the Modeltest code.");
    fprintf(fp, "\n______________________________________________________________________\n\n");
}

/********************* PrintDate ***********************/
static void PrintDate (FILE *fp)
{
    time_t now;
    char *date;

    now = time(NULL);
    date = ctime(&now);
    fprintf(fp, "%s",date);
}

/************** CheckNA **********************/
/*
    If value is NA prints "-"

*/
//static char *CheckNA (double value)
//{
//    char *string;
//
//    string = (char*) calloc (100, sizeof (char));
//    if (value == NA) {
//        return "  -  ";
//    }
//    else {
//        sprintf (string, "%8.4f", value);
//        return string;
//    }
//}

/*********************** PrintUsage ********************************/
static void PrintUsage()
{
    fprintf(stderr,"\n HELP \n");
    fprintf(stderr,"\nMrModeltest is a program for comparing models of evolution using likelihood ");
    fprintf(stderr,"ratio tests and the AIC criterion. The input are log likelihood scores. ");
    fprintf(stderr,"You can input raw scores or a Paup matrix resulting from the execution ");
    fprintf(stderr,"of the provided block of Paup commands (MrModelblock).");
    fprintf(stderr,"\nThe program can also enter in a calculator mode for obtaining ");
    fprintf(stderr,"the P-value associated with the log likelihood ratio statistic for two given ");
    fprintf(stderr,"scores or the AIC value for entered scores.\n");
    fprintf(stderr,"\nNOTE --- MrModeltest only tests 24 models (Modeltest uses 56). However ");
    fprintf(stderr,"all of the 24 models can be specified in MrBayes (version 3).");
    fprintf(stderr," MrModeltest also use (by default) four different hierarchies when conducting the");
    fprintf(stderr," likelihood ratio tests. The hierarchies are described in detail by");
    fprintf(stderr," Posada & Crandall. 2001. Systematic Biology, 50:580-601 (Figure 4).");
    fprintf(stderr,"\n\nJC:    Jukes and Cantor 1969");
    fprintf(stderr,"\nK80:   Kimura 2-parameters, Kimura 1980 (also known as K2P)");
    fprintf(stderr,"\nSYM:   Symmetrical model, Zharkikh 1994");
    fprintf(stderr,"\nF81:   Felsenstein 1981");
    fprintf(stderr,"\nHKY:   Hasegawa-Kishino-Yano 1985");
    fprintf(stderr,"\nGTR:   General time reversible, Rodriguez et al 1990 (also known as REV)");
    fprintf(stderr,"\nI:     invariable sites");
    fprintf(stderr,"\nG:     gamma distribution");
    fprintf(stderr,"\n\nUsage:");
    fprintf(stderr,"\n         -? : help");
    fprintf(stderr,"\n         -2 : use alternative hLRT2 hierarchy (starting with GTR+I+G vs. SYM+I+G)");
    fprintf(stderr,"\n         -3 : use alternative hLRT3 hierarchy (starting with JC vs. JC+G)");
    fprintf(stderr,"\n         -4 : use alternative hLRT4 hierarchy (starting with GTRIG vs. GTRI)");
    fprintf(stderr,"\n         -a : alpha level (e.g. -a0.01)");
    fprintf(stderr,"\n         -d : debug level (e.g. -d2)");
    fprintf(stderr,"\n         -f : input from a file for obtaining AIC values");
    fprintf(stderr,"\n         -h : help");
    fprintf(stderr,"\n         -i : AIC calculator mode");
    fprintf(stderr,"\n         -l : LRT calculator mode");
    fprintf(stderr,"\n         -n : sample size or number of characters (all or just variable). Forces the use of AICc");
    fprintf(stderr,"\n         -t : number of taxa. Forces to include branch lengths as parameters");
    fprintf(stderr,"\n         -v : prints version number");
    /*fprintf(stderr,"\n         -w : confidence interval for averaging (e.g., -w0.95) (default is w=1.0)");*/
    fprintf(stderr,"\n\nUNIX/MACOSX/WIN usage: mrmodeltest2 [-d -a -c -t -2 -3 -4 -l -i -f -? -h] < mrmodel.scores > outfile\n\n");
    /*fprintf(stderr,"\n\nHit return to close this window!\n\n");*/    /*For windows.*/
}

/********************** Allocate *************************/
static void Allocate()
{
    modelhLRT = (char*) calloc (10, sizeof (char));
    modelhLRT2 = (char*) calloc (10, sizeof (char));
    modelhLRT3 = (char*) calloc (10, sizeof (char));
    modelhLRT4 = (char*) calloc (10, sizeof (char));
    modelAIC = (char*) calloc (10, sizeof (char));
    model = (ModelSt*) calloc (NUM_MODELS, sizeof (ModelSt));
    order = (ModelSt*) calloc (NUM_MODELS, sizeof (ModelSt));
}

/*********************** Free ***************************/
static void Free()
{
    free (modelhLRT);
    free (modelhLRT2);
    free (modelhLRT3);
    free (modelhLRT4);
    free (modelAIC);
    free (model);
    free (order);
}

