#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct pixel
{
    unsigned char R;
    unsigned char G;
    unsigned char B;
};

struct correlation
{
    int i;
    int j;
    int sablon;
    double value;

};

unsigned int * xorshift32 (unsigned int seed, unsigned int width, unsigned int height)
{
    unsigned int * arr_of_R =(unsigned int *) malloc(sizeof(int)*2*width*height);

    for (int i=1; i<width*height*2; i++)
    {
        seed = seed ^ seed << 13;
        seed = seed ^ seed >> 17;
        seed = seed ^ seed << 5;
        arr_of_R[i] = seed;

    }
    return arr_of_R;

}

struct pixel * linearization_imagine (char *image_bmp)
{
    unsigned int width, height;

    FILE *f = fopen (image_bmp, "rb");

    if (f==NULL)
    {
        printf("Nu exista fisierul pentru criptare");
        exit(0);
        return NULL;
    }

    int padding;

    fseek(f,18,SEEK_SET);
    fread(&width, sizeof(int), 1, f);
    fseek(f,22,SEEK_SET);
    fread(&height, sizeof (int), 1,f);


    if (3*width % 4 == 0)
        padding = 0;
    else
        padding = 4 - (3*width) %4;

    printf("padding = %d\n", padding);

    int k;

    printf("latime = %u ", width);
    printf("\ninaltime = %u \n", height);

    struct pixel * l_i =(struct pixel *) malloc (width*height*sizeof(struct pixel));
    unsigned int n=0;
    unsigned char p[3];
    unsigned int i;
    fseek(f,54,SEEK_SET);
    while (fread(&p,3,1,f)==1)
    {
        l_i[n].R = p[2];
        l_i[n].G = p[1];
        l_i[n].B = p[0];

        if ((n%width) == (width - 1) && n!=0)
        {
            fseek(f,padding, SEEK_CUR);
        }
        n++;

    }

    switching_pixels(l_i, width,height);
    fclose(f);
    return l_i;
}

void switching_pixels (struct pixel * linear_image,  unsigned int width, unsigned int height)
{
    unsigned int n = width * height;
    int i, j;
    struct pixel aux;
    for (i=0; i<=n/2; i++)
    {
        aux.R = linear_image[n-i].R;
        aux.G = linear_image[n-i].G;
        aux.B = linear_image[n-i].B;

        linear_image[n-i].R = linear_image[i].R;
        linear_image[n-i].G = linear_image[i].G;
        linear_image[n-i].B = linear_image[i].B;

        linear_image[i].R = aux.R;
        linear_image[i].G = aux.G;
        linear_image[i].B = aux.B;

    }

    for (i=0; i<height; i++)
    {
        for (j=0; j<=width/2; j++)
        {
            aux.R = linear_image[i*width + j].R;
            aux.G = linear_image[i*width + j].G;
            aux.B = linear_image[i*width + j].B;

            linear_image[i*width + j].R = linear_image[i*width + width - j].R;
            linear_image[i*width + j].G = linear_image[i*width + width - j].G;
            linear_image[i*width + j].B = linear_image[i*width + width - j].B;

            linear_image[i*width + width - j].R = aux.R;
            linear_image[i*width + width - j].G = aux.G;
            linear_image[i*width + width - j].B = aux.B;
        }
    }
}

unsigned int * random_permutation (unsigned int width, unsigned int height, unsigned int * R)
{
    unsigned int * arr = malloc (height*width*sizeof(int));
    int i;

    for (i = 0; i < height*width; i++)
    {
        arr[i]=i;
    }

    int j;

    for (i = height*width-1; i>0; i--)
    {
        j = R[height*width - i]%(i+1);
        int aux = arr[i];
        arr[i] = arr[j];
        arr[j] = aux;
    }

    return arr;

}

void pixel_permutation(struct pixel * linear_image, unsigned int *arr_of_permutations, struct pixel *li, unsigned int width, unsigned int height)
{
     int i;
    for (i = 0; i < width*height; i++)
    {
        li[i].R = linear_image[i].R;
        li[i].G = linear_image[i].G;
        li[i].B = linear_image[i].B;
    }

    for (i = 0; i < width*height; i++)
    {
        linear_image[arr_of_permutations[i]].R = li[i].R;   ///Red
        linear_image[arr_of_permutations[i]].G = li[i].G;   ///Green
        linear_image[arr_of_permutations[i]].B = li[i].B;   ///Blue
    }

}

struct pixel *  change (struct pixel * linear_image, unsigned int * arr_of_R,  char* orignal_bmp, char * imagine_bmp, char * fisiertxt)
{
    FILE *f = fopen(orignal_bmp, "rb");
    FILE *fin = fopen(fisiertxt, "r");
    if (f == NULL || fin == NULL)
    {
        printf("nu exista fisierul original");
        exit(0);
    }
    unsigned int SV;
    unsigned int width;
    unsigned int height;

    fseek(f,18,SEEK_SET);
    fread(&width, sizeof(int), 1, f);
    fseek(f,22,SEEK_SET);
    fread(&height, sizeof (int), 1,f);


    fscanf(fin,"%u",&SV);
    fscanf(fin,"%u",&SV);
    printf("%u", SV);

    unsigned char *p;
    unsigned char *p2;
    p=&SV;
    struct pixel * C = (struct pixel *) malloc (width*height*sizeof(struct pixel));
    int i;
    for (i = 0; i < width*height; i++)
    {
        p2 = &arr_of_R[width*height + i];

        if (i == 0)
        {
            C[i].B = linear_image[i].B ^(*p)^(*p2);     ///Blue
            p = p + 1;
            p2 = p2+1;
            C[i].G = linear_image[i].G ^(*p)^(*p2);     ///Green
            p = p + 1;
            p2 = p2+1;
            C[i].R = linear_image[i].R ^(*p)^(*p2);     ///Red
        }
        else
        {
            C[i].B = (linear_image[i].B) ^(C[i-1].B)^(*p2);      ///Blue
            p2 = p2+1;
            C[i].G = (linear_image[i].G) ^(C[i-1].G)^(*p2);      ///Green
            p2 = p2+1;
            C[i].R = (linear_image[i].R) ^(C[i-1].R)^(*p2);      ///Red
        }

    }
    fclose(fin);
    fclose(f);
    afis(C,width,height,orignal_bmp,imagine_bmp);
    return C;
}

void rechange(struct pixel *D, struct pixel *C, unsigned int *arr_of_R,unsigned int SV, unsigned int width, unsigned int height)
{
    unsigned char *p = &SV;
    unsigned char *p2;
    for (int i = 0; i < width*height; i++)
    {

        unsigned char *p2;
        p2 = &arr_of_R[width*height + i];

        if (i == 0)
        {
            D[i].B = C[i].B ^(*p)^(*p2);     ///Blue
            p = p + 1;
            p2 = p2+1;
            D[i].G = C[i].G ^(*p)^(*p2);     ///Green
            p = p + 1;
            p2 = p2+1;
            D[i].R = C[i].R ^(*p)^(*p2);     ///Red
        }
        else
        {
            D[i].B = (C[i].B) ^(C[i-1].B)^(*p2);      ///Blue
            p2 = p2+1;
            D[i].G = (C[i].G) ^(C[i-1].G)^(*p2);      ///Green
            p2 = p2+1;
            D[i].R = (C[i].R) ^(C[i-1].R)^(*p2);      ///Red
        }

    }
}

void chi_test (struct pixel* linear_image, unsigned int width, unsigned int height)
{
    unsigned int *red_channel;
    unsigned int *green_channel;
    unsigned int *blue_channel;

    ///Creates space in HEAP for the arrays for the chi test
    red_channel = (unsigned int*)malloc(sizeof(unsigned int)*256);
    green_channel = (unsigned int*)malloc(sizeof(unsigned int)*256);
    blue_channel = (unsigned int*)malloc(sizeof(unsigned int)*256);

    double red = 0, blue = 0, green = 0, formula_value;
    ///Sets the 3 arrays of the channels to 0
    for (int i=0; i<256; i++)
    {
        red_channel[i] = green_channel[i] = blue_channel[i] = 0;
    }

    ///Puts the frequency of the values of the channels in the 3 specified arrays
    for (int i=0; i<width*height; i++)
    {
        red_channel[linear_image[i].R]++;
        green_channel[linear_image[i].G]++;
        blue_channel[linear_image[i].B]++;
    }
    ///Sets the values of the channels with that straneg formula
    formula_value = (double)(width*height/256);

    ///The chi-test (before cyphering)
    for (int i=0; i<256; i++)
    {
        double aux = (double)((double)((double)(red_channel[i] - formula_value)*(double)(red_channel[i] - formula_value))/formula_value);
        red += aux;

        aux = (double)((double)((double)(green_channel[i] - formula_value)*(double)(green_channel[i] - formula_value))/formula_value);
        green += aux;

        aux = (double)((double)((double)(blue_channel[i] - formula_value)*(double)(blue_channel[i] - formula_value))/formula_value);
        blue += aux;
    }
    ///Showing the results
    printf("rosu: %f\nverde: %f\nalbastru: %f\n", red,green,blue);
    free(red_channel);
    free(blue_channel);
    free(green_channel);
}

void afis (struct pixel * D, unsigned int width, unsigned int height, char * original_bmp,char * imagine_bmp)
{
    FILE *fout = fopen (imagine_bmp, "wb+");
    FILE *f = fopen (original_bmp, "rb");

    if (f==NULL)
    {
        printf("Nu exista fisierul pentru criptare");
        exit(0);
        return NULL;
    }
    unsigned char *buffer = (unsigned char *)malloc(sizeof(unsigned char)*54);                   ///Array of char for containing the header of the BMP

    ///Copies the header from the original picture into the second one
    fseek(f,0,SEEK_SET);
    fread(buffer,54,1,f);
    fwrite(buffer,54,1,fout);
    fclose(f);

    unsigned int n = 0;
    unsigned char d = 0;
    int padding;

    if (width % 4 == 0)
        padding = 0;
    else
        padding = 4 - (3*width) %4;

    while (n < width*height)
    {
        fwrite(&D[n].B, 1,1,fout);
        fwrite(&D[n].G, 1,1,fout);
        fwrite(&D[n].R, 1,1,fout);

        if ((n%width) ==(width -1) && n!=0)
        {
            int k=0;
            unsigned char d = '0';
            if (k!=padding)
                while (k<padding)
                {
                    fwrite(&d, 1, 1, fout);
                    k++;
                }
        }
        n++;
    }
    free(buffer);
    fclose(fout);
}

void grayscale(struct pixel * linear_image, unsigned int width, unsigned int height)
{
    unsigned char aux;
    for (int i = 0; i < width*height; i++)
    {
        aux = 0.299*linear_image[i].R + 0.587*linear_image[i].G + 0.114*linear_image[i].B;
        linear_image[i].R = linear_image[i].G = linear_image[i].B = aux;
    }
}

void color(struct pixel *D, unsigned int x, unsigned int y, unsigned int width, unsigned int height, struct pixel RGB)
{
    for (int i = y; i<y+11; i++)
    {
        D[x*width + i].R = D[(x+15)*width + i].R = RGB.R;
        D[x*width + i].G = D[(x+15)*width + i].G = RGB.G;
        D[x*width + i].B = D[(x+15)*width + i].B = RGB.B;
    }

    for (int i = 0; i<15; i++)
    {
        D[(x+i)*width + y].R = D[(x+i)*width + y + 11].R = RGB.R;
        D[(x+i)*width + y].G = D[(x+i)*width + y + 11].G = RGB.G;
        D[(x+i)*width + y].B = D[(x+i)*width + y + 11].B = RGB.B;
    }
}

double temp_plate (struct pixel *D, struct pixel * t, int i, int j,double ps,unsigned int width, unsigned int height, unsigned int template_width, unsigned int template_height)
{
    double medie_template = 0;
    double medie_image = 0;

    for (int z = 0; z<template_height*template_width; z++)
    {
        medie_template += t[z].R;
    }
    medie_template = medie_template / (float)(template_height*template_width);

//printf("medie template = %f\n", medie_template);


    double sigma_template = 0;
    double sigma_image = 0;
    double aux_template =0;
    double aux_image =0;

    for (int z = 0; z<15; z++)
    {
        for (int q = 0; q< 11; q++)
            medie_image += D[(z+i)*width + q + j].R;
    }

    medie_image = medie_image / (float)(165);
    //printf("medie image = %f\n", medie_image);
    aux_template =0;
    aux_image =0;
    for (int x = 0; x < template_height; x++)
    {
        for (int y = 0; y < template_width; y++)
        {
            aux_template +=(double)(t[x*template_width+y].R - medie_template)*(double)(t[x*template_width+y].R - medie_template);
            aux_image +=(double)(D[(i+x)*width+y+j].R - medie_image)*(double)(D[(i+x)*width+y+j].R - medie_image);
        }
    }
    aux_template = aux_template / (double)(template_height*template_width-1);
    sigma_template = sqrt(aux_template);


    aux_image = aux_image / (double)(template_height*template_width-1);
    sigma_image = sqrt(aux_image);

    aux_image = 0;
    aux_template = 0;
    double aux_aux = 0;
    for (int x = 0; x < template_height; x++)
    {
        for (int y = 0; y < template_width; y++)
        {
            aux_template =(double)(t[x*template_width+y].R - medie_template);
            aux_image =(double)(D[(i+x)*width+y+j].R - medie_image);
            aux_aux+=(aux_image*aux_template)/(sigma_image*sigma_template);
        }
    }
    aux_aux = aux_aux/(double)(template_height*template_width);
    if (aux_aux>=ps)
    {
        printf("%d %d %f  ", i, j, aux_aux);
         return aux_aux;

    }


   return -1;

}

int compare(const void *a, const void *b)
{
    struct correlation va = *(struct correlation *)a;
    struct correlation vb = *(struct correlation *)b;
    if(va.value<vb.value)
        return 1;
    if(va.value>vb.value)
        return -1;
    return 0;
}

double suprapunere (struct correlation di, struct correlation dj, int ariei, int ariej)
{

    if (abs(di.i - dj.i)>15 || abs(di.j - dj.j)>11)
        return 0;
    printf("%d %d : %d %d ", di.i, di.j, dj.i, dj.j);

    double arie_i, arie_j, arie_ij;
    arie_i = ariei;
    arie_j = ariej;

    int min_x, max_x, min_y, max_y;

    if (di.i > dj.i)
    {

        min_x = dj.i;
        max_x = di.i;
    }
    else
    {
        min_x = di.i;
        max_x = dj.i;
    }
    if (di.j > dj.j)
    {

        min_y = dj.j;
        max_y = di.j;
    }
    else
    {
        min_y = di.j;
        max_y = dj.j;
    }

    if (min_x == max_x && min_y == max_y)
        arie_ij = 165;
    else if (min_x == max_x)
    {
        arie_ij = 15*(11 - (max_y - min_y));
    }
    else if (min_y == max_y)
        arie_ij = (15 - ( max_x - min_x))*11;
    else
        arie_ij = (15  - (max_x - min_x))*(11 - (max_y - min_y));
    double suprapuneri = arie_ij/(arie_i + arie_i - arie_ij);
    printf("%lf \n", suprapuneri);
    return suprapuneri;
}

void remove_correlation (struct correlation *D_correlation, int m, int ariei, int ariej)
{

    int i = 0;
    int j = 0;
    while (i<m)
    {
        j=i+1;
        while (j<m)
        {
            if ( suprapunere(D_correlation[i],D_correlation[j], ariei, ariej)>0.2)
            {
                for (int z=j; z<m-1; z++)
                    D_correlation[z] = D_correlation[z+1];
                m--;
                j--;
            }
            j++;
        }
        i++;
    }
}

int main()
{
    char *original_bmp = "test.bmp";                                   ///Original bmp file for cyphering
    char *imagine_bmp = "imagine_criptata_si_decriptata.bmp";          ///Destination bmp file after cyphering and decrypting
    char *template_bmp = "imagine_template.bmp";                        ///Destination bmp file after template
    unsigned int height, width;                                         ///Variables taken from the buffer for the dimensions of the BMP
    unsigned int *arr_of_R;                                             ///Array with random generated numbers from the R0 seed
    unsigned int *arr_of_permutations;                                  ///Array with random permutations for the pixels in the image
    unsigned int R0 = 123456789;                                        ///First Seed for XORSHIFT32
    unsigned int SV = 987654321;                                        ///Secret Value
    struct pixel *linear_image;                                         ///Array of holding the pixels of the picture in one single array
    struct pixel *li;                                                   ///Copy of the linear_image variable (the pixels in this variable won't be changed)
    struct pixel *C;                                                    ///Array that holds the cyphered array of pixels
    struct pixel *D;                                                    ///Array that hold the decrypted array of pixels
    unsigned char * p = &SV;                                            ///Variable that stores the bytes of the SV
    unsigned char b[3];                                                 ///Variable that holds a pixel (BGR)
    unsigned int *arr_of_cpermutations;                                 ///Array the holds the revers permutation
    unsigned char * p2;
    unsigned int n=0;

    ///Finds the dimensions of the BMP picture in the header
    FILE *f = fopen (original_bmp, "rb");
    fseek(f,18,SEEK_SET);
    fread(&width, sizeof(int), 1, f);
    fseek(f,22,SEEK_SET);
    fread(&height, sizeof (int), 1,f);
    printf("latime = %u ", width);
    printf("\ninaltime = %u \n", height);
    fclose(f);



    ///Copies the address of the dynamic created array of random generated numbers using the R0 seed in this variable
    arr_of_R = xorshift32(R0, width, height);

    ///Copies the address of the dynamic created array of random permutations for the pixels
    arr_of_permutations = random_permutation(width, height, arr_of_R);

    ///Copies the address of the dynamic created array of the pixels of the image
    linear_image = linearization_imagine(original_bmp);

    li = (struct pixel *) malloc (width*height*sizeof(struct pixel));

    struct pixel * D_original = (struct pixel *) malloc (width*height*sizeof(struct pixel));

    arr_of_cpermutations = (unsigned int*)malloc(sizeof(unsigned int)*width*height);

    D = (struct pixel *) malloc (width*height*sizeof(struct pixel));

    for (int i = 0; i<width*height; i++)
        D_original[i] = linear_image[i];

    ///Shows the result of the chi-test before cyphering the image
    printf("\nInainte de criptare\n");
    chi_test (linear_image, width, height);

    ///Interchanges the pixels in the image
    pixel_permutation(linear_image, arr_of_permutations, li, width,height);

    ///Creates the reversed permutation array
    for (int i=0; i<width*height; i++)
    {
        arr_of_cpermutations[arr_of_permutations[i]] = i;
    }

    ///Cyphers the pixels
    C = change(linear_image,arr_of_R,original_bmp,"criptat.bmp","secret_key.txt");

    ///Chi-test after cyphering
    printf("\nDupa criptare\n");
    chi_test (C, width, height);

    ///Decrypts the pixels
    rechange(D,C,arr_of_R,SV,width,height);

    ///Re-interchanges the pixels in the image
    pixel_permutation(D, arr_of_cpermutations, li, width,height);




///Here ends the cyphering part

///Here starts the template part


    ///Reputs the pixels just like before the linearization
    switching_pixels(D,width,height);

    ///Puts the pixels in the new .bmp folder
    afis(D,width, height,original_bmp,imagine_bmp);

    ///Puts the pixel just like in the linearization
    switching_pixels(D,width,height);

    ///Gray-scales the bmp file
    grayscale(D,width,height);

    char zero_binary[] = "cifra0.bmp";
    char one_binary[] = "cifra1.bmp";
    char two_binary[] = "cifra2.bmp";
    char three_binary[] = "cifra3.bmp";
    char four_binary[] = "cifra4.bmp";
    char five_binary[] = "cifra5.bmp";
    char six_binary[] = "cifra6.bmp";
    char seven_binary[] = "cifra7.bmp";
    char eight_binary[] = "cifra8.bmp";
    char nine_binary[] = "cifra9.bmp";

    ///The arrays that will hold the linearization of the template bmp files
    struct pixel * zero_bmp = linearization_imagine(zero_binary);
    struct pixel * one_bmp = linearization_imagine(one_binary);
    struct pixel * two_bmp = linearization_imagine(two_binary);
    struct pixel * three_bmp = linearization_imagine(three_binary);
    struct pixel * four_bmp = linearization_imagine(four_binary);
    struct pixel * five_bmp = linearization_imagine(five_binary);
    struct pixel * six_bmp = linearization_imagine(six_binary);
    struct pixel * seven_bmp = linearization_imagine(seven_binary);
    struct pixel * eight_bmp = linearization_imagine(eight_binary);
    struct pixel * nine_bmp = linearization_imagine(nine_binary);

    ///Gray-scales the templates
    grayscale(zero_bmp,11,15);
    grayscale(one_bmp,11,15);
    grayscale(two_bmp,11,15);
    grayscale(three_bmp,11,15);
    grayscale(four_bmp,11,15);
    grayscale(five_bmp,11,15);
    grayscale(six_bmp,11,15);
    grayscale(seven_bmp,11,15);
    grayscale(eight_bmp,11,15);
    grayscale(nine_bmp,11,15);

    ///creates an array that holds the colors for the windows
    struct pixel colors[10];
    colors[0].R = 255;
    colors[0].G = 0;
    colors[0].B = 0;

    colors[1].R = 255;
    colors[1].G = 255;
    colors[1].B = 0;

    colors[2].R = 0;
    colors[2].G = 255;
    colors[2].B = 0;


    colors[3].R = 0;
    colors[3].G = 255;
    colors[3].B = 255;

    colors[4].R = 255;
    colors[4].G = 0;
    colors[4].B = 255;

    colors[5].R = 0;
    colors[5].G = 0;
    colors[5].B = 255;

    colors[6].R = 192;
    colors[6].G = 192;
    colors[6].B = 192;

    colors[7].R = 255;
    colors[7].G = 140;
    colors[7].B = 0;

    colors[8].R = 128;
    colors[8].G = 0;
    colors[8].B = 128;

    colors[9].R = 128;
    colors[9].G = 0;
    colors[9].B = 0;

    ///Dimensions of the templates
    unsigned int template_height, template_width;
    f = fopen (zero_binary,"rb");
    if (f == NULL)
    {
        return 1;
    }
    struct correlation * D_correlation = (struct correlation *)malloc(sizeof(struct correlation)*width*height*10);
    fseek(f,18,SEEK_SET);
    fread(&template_width, sizeof(int), 1, f);
    fseek(f,22,SEEK_SET);
    fread(&template_height, sizeof (int), 1,f);
    fclose(f);
    long int m = 0;
    double ps = 0.5;

    double *max = (double *) malloc(sizeof(double)*10);

    ///Goes through the image and puts in the correlation array all the position of the pixels, the correlation value and the color that must be used
    for (int i = 0; i <= height - template_height; i++)
    {
        for (int j = 0; j <= width - template_width; j++)
        {
            max[0] = temp_plate(D,zero_bmp,i,j,ps,width,height,template_width,template_height);
            max[1] = temp_plate(D,one_bmp,i,j,ps,width,height,template_width,template_height);
            max[2] = temp_plate(D,two_bmp,i,j,ps,width,height,template_width,template_height);
            max[3] = temp_plate(D,three_bmp,i,j,ps,width,height,template_width,template_height);
            max[4] = temp_plate(D,four_bmp,i,j,ps,width,height,template_width,template_height);
            max[5] = temp_plate(D,five_bmp,i,j,ps,width,height,template_width,template_height);
            max[6] = temp_plate(D,six_bmp,i,j,ps,width,height,template_width,template_height);
            max[7] = temp_plate(D,seven_bmp,i,j,ps,width,height,template_width,template_height);
            max[8] = temp_plate(D,eight_bmp,i,j,ps,width,height,template_width,template_height);
            max[9] = temp_plate(D,nine_bmp,i,j,ps,width,height,template_width,template_height);

            D_correlation[m].value = max[0];
            D_correlation[m].sablon = 0;
            D_correlation[m].i = i;
            D_correlation[m].j = j;
            m++;
            D_correlation[m].value = max[1];
            D_correlation[m].sablon = 1;
            D_correlation[m].i = i;
            D_correlation[m].j = j;
            m++;
            D_correlation[m].value = max[2];
            D_correlation[m].sablon = 2;
            D_correlation[m].i = i;
            D_correlation[m].j = j;
            m++;
            D_correlation[m].value = max[3];
            D_correlation[m].sablon = 3;
            D_correlation[m].i = i;
            D_correlation[m].j = j;
            m++;
            D_correlation[m].value = max[4];
            D_correlation[m].sablon = 4;
            D_correlation[m].i = i;
            D_correlation[m].j = j;
            m++;

            D_correlation[m].value = max[5];
            D_correlation[m].sablon = 5;
            D_correlation[m].i = i;
            D_correlation[m].j = j;
            m++;

            D_correlation[m].value = max[6];
            D_correlation[m].sablon = 6;
            D_correlation[m].i = i;
            D_correlation[m].j = j;
            m++;

            D_correlation[m].value = max[7];
            D_correlation[m].sablon = 7;
            D_correlation[m].i = i;
            D_correlation[m].j = j;
            m++;
            D_correlation[m].value = max[8];
            D_correlation[m].sablon = 8;
            D_correlation[m].i = i;
            D_correlation[m].j = j;
            m++;
            D_correlation[m].value = max[9];
            D_correlation[m].sablon = 9;
            D_correlation[m].i = i;
            D_correlation[m].j = j;
            m++;
        }

    }
    printf ("m = %d\n", m);
    ///Sorts the correlation array
    qsort(D_correlation,m,sizeof(struct correlation), compare);

    ///Removes the correlations that are smaller than ps
    for (int z=0; z<m; z++)
    {
        if (D_correlation[z].value<ps)
            m=z;
    }
    printf ("m = %d\n", m);

    ///Goes through the correlation array and removes the objects that have the correlation value smaller than the current correlation and their surfaces are touching each other
    remove_correlation(D_correlation,m, template_height*template_width, template_height*template_width);

    for (int i = 0; i<m; i++)
    {
        color(D_original,D_correlation[i].i,D_correlation[i].j,width,height,colors[D_correlation[i].sablon]);
    }


    ///Reputs the pixels just like before the linearization
    switching_pixels(D_original,width,height);

    ///Puts the pixels in the new .bmp folder
    afis(D_original,width, height,original_bmp,template_bmp);

    free(arr_of_cpermutations);
    free(arr_of_permutations);
    free(arr_of_R);
    free(linear_image);
    free(D);
    free(D_correlation);
    free(max);
    free(D_original);
    free(li);
    free(C);
    return 0;
}
