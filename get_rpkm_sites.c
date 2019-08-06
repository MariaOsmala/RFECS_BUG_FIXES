/************************************************************************************/

#include "stdio.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "float.h"

#define MAX_CHROM_NUM 24       //maximum chromosome number

typedef struct
{
    char chromName[100];
    int chromSize;
}CHROM_INFO;

//Structure defining the histone modification regions

CHROM_INFO chromInfo[MAX_CHROM_NUM];

//total number of bins
int totalBinNum;

int chromNum;

int ChromToIndex(char *chrom);
char *IndexToChrom(int index, char *chrom, int len);
int GetChromInfo(char *fileName);
int GetRPKM(char *siteFileName, char *RPKMFilename,int width, int bins,char *projectName);

int main(int argc, char* argv[])
{
    char RPKMFileName[1000],siteFileName[1000],projectName[1000],chromFileName[1000];
    int width,bins;
    int count=0;
    if (argc!=7)
    {
        printf("Usage: <width> <number of bins> <sites> <RPKM file> <chromsome coordinate file><project name>\n");
        return 0;
    }

    width=atoi(argv[1]);
    bins=atoi(argv[2]);
    strcpy(siteFileName,argv[3]);
    strcpy(RPKMFileName,argv[4]);
    strcpy(chromFileName,argv[5]);
    strcpy(projectName,argv[6]);


    //Read chromosome description file
    if (!GetChromInfo(chromFileName))
    {
        printf("chromosome description file is not valid\n");
        return 0;
    }

    //Read tag files of L1 and L2, and define the histone modification regions
    GetRPKM(siteFileName,RPKMFileName,width,bins,projectName);



    return 0;
}

//transform chromosome name the chromosome index
//
//
//
//
//
int ChromToIndex(char *chrom)
{
    int i;

    for (i=0;i<chromNum;i++)
    {
        if (!strcmp(chrom, chromInfo[i].chromName))
        {
            break;
        }
    }

    if (i<chromNum)
    {
        return i;
    }
    else
    {
        return -1;
    }
}

//transfrom chromosome index to chromosome name
char *IndexToChrom(int index, char *chrom, int len)
{
    if ((index<0)||(index>=chromNum))
    {
        return 0;
    }

    if (strlen(chromInfo[index].chromName)>=len)
    {
        return 0;
    }

    strcpy(chrom, chromInfo[index].chromName);

    return chrom;
}

//Read the chromosome description file
int GetChromInfo(char *fileName)
{
    FILE *fh;
    char tmpStr[1000];

    fh = (FILE *)fopen(fileName, "r");

    if (!fh)
    {
        return -1;
    }

    chromNum = 0;

    fscanf(fh, "%s", tmpStr);

    while (!feof(fh))
    {
        strcpy(chromInfo[chromNum].chromName, tmpStr);
        fscanf(fh, "%s", tmpStr);
        chromInfo[chromNum].chromSize = atoi(tmpStr);
        fscanf(fh, "%s", tmpStr);
        chromNum++;
    }

    fclose(fh);
    return chromNum;
}

//Read the tag files of L1 and L2, and determine the histone modification sites
int GetRPKM(char *siteFileName, char *RPKMFilename,int width, int bins,char *projectName)
{
    FILE *fh,*fh2;
    char tmpStr[1000];
    char tmpChrom[100];
    int tmpPos,tmpIndex,tmpChrom2;
    int tmpStart, tmpEnd;
    char tmpStrand;
    int i,j,k;
    float tmpVal;

    int BIN_SIZE=(2*width)/bins;
    //binCounts and binCounts2 store the fragment counts in each bin. mask=1 flags histone modification site

    float **binVals;
    binVals = (float **)malloc(chromNum*sizeof(int *));

    //Initialize the arrays
    totalBinNum = 0;

    for (i=0;i<chromNum;i++)
    {
        totalBinNum += chromInfo[i].chromSize/BIN_SIZE+1;
        binVals[i] = (float *)malloc((chromInfo[i].chromSize/BIN_SIZE+1)*sizeof(float));
        memset(binVals[i], 0, (chromInfo[i].chromSize/BIN_SIZE+1)*sizeof(float));
    }

    fh = (FILE *)fopen(RPKMFilename,"r");

    if (!fh)
    {
        return -1;
    }


    fscanf(fh, "%d",&tmpChrom2);
    tmpIndex=tmpChrom2-1;
    while (!feof(fh))
    {

        if (tmpIndex<0)
        {
            fscanf(fh, "%d", &tmpPos); fscanf(fh, "%f", &tmpVal);
            fscanf(fh, "%d", &tmpChrom2);
            tmpIndex=tmpChrom2-1;
            continue;
        }

        fscanf(fh, "%d", &tmpPos);fscanf(fh, "%f", &tmpVal);
        if ((tmpPos>=chromInfo[tmpIndex].chromSize)||(tmpPos<0))
        {

            fscanf(fh, "%d", &tmpChrom2);
            tmpIndex=tmpChrom2-1;
            continue;
        }
        binVals[tmpIndex][tmpPos/BIN_SIZE]=tmpVal;
        fscanf(fh,"%d", &tmpChrom2);tmpIndex=tmpChrom2-1;
    }

    printf("Done reading RPKM file\n");

    fclose(fh);
    char binFileName[1000];
    sprintf(binFileName,"%s.list",projectName);

    fh2 = (FILE *)fopen(binFileName,"w");

    fh = (FILE *)fopen(siteFileName, "r");
    if (!fh)
    {
        return -1;
    }



    fscanf(fh, "%s", tmpChrom);

    while (!feof(fh))
    {
        tmpIndex = ChromToIndex(tmpChrom);
        if (tmpIndex<0)
        {
            fscanf(fh, "%d", &tmpPos); fscanf(fh, "%s", &tmpStr);
            fscanf(fh, "%s", tmpChrom);
            continue;
        }

        fscanf(fh, "%d", &tmpPos);fscanf(fh, "%s", &tmpStr);
        tmpStrand = tmpStr[0];
        int abc = tmpPos/BIN_SIZE;
        if(tmpStrand=='+')
        {

            tmpStart=(abc*BIN_SIZE)-width;
            tmpEnd=(abc*BIN_SIZE)+width-BIN_SIZE;
        }
        if(tmpStrand=='-')
        {
            tmpStart=(abc*BIN_SIZE)-width+BIN_SIZE;
            tmpEnd=(abc*BIN_SIZE)+width;
        }
        if ((tmpPos>=chromInfo[tmpIndex].chromSize)||(tmpPos<0))
        {
            fscanf(fh, "%s", tmpChrom);
            continue;
        }

        fscanf(fh,"%s", tmpChrom);
        for(i=tmpStart;i<tmpEnd;i=i+BIN_SIZE)
        {
            fprintf(fh2,"%f\t",binVals[tmpIndex][i/BIN_SIZE]);
        }
        fprintf(fh2,"%f\n",binVals[tmpIndex][tmpEnd/BIN_SIZE]);
    }



    return 0;
}




