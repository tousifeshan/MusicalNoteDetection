//
// Created by Tousif Ahmed on 2/10/16.
//
#include <iostream>
using namespace std;

typedef enum {NOTEHEAD=0, QUARTERREST=1, EIGHTHREST=2} Type;
class DetectedSymbol {
public:
    int row, col, width, height;
    Type type;
    char pitch;
    double confidence;
};
vector<int> staffs;
int avg_space=11;

int copy_edge(int val, int boundary);
char get_pitch(int _top, int _left, int _bottom, int _right);
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
    SDoublePlane output(input.rows(), input.cols());


    //Creates a new image with same size as the input initialized to all 0s (black)
    SDoublePlane temp_output(input.rows(), input.cols());


    int x, y, fx, fy, input_x, input_y;
    double accum;

    for(x = 0; x < input.rows(); x++) {
        for(y = 0; y < input.cols(); y++) {
            accum = 0.0;
            for(fy = 0; fy < col_filter.rows(); fy++) {
                //input_y = y + (fy - col_filter.cols() / 2);
                input_y = copy_edge(y + (fy - col_filter.rows() / 2), input.cols());
                accum += col_filter[fy][0] * input[x][input_y];
            } // last filter for-loop
            temp_output[x][y] = accum;
        }
    }
    // apply the hoirz column to the image up then down, left to right
    for(y = 0; y < temp_output.cols(); y++) {
        for(x = 0; x < temp_output.rows(); x++) {
            accum = 0.0;
            for(fx = 0; fx < row_filter.cols(); fx++) {
                //input_x = x + (fx - row_filter.rows()/2);
                input_x = copy_edge(x + (fx - row_filter.cols()/2), temp_output.rows());
                accum += row_filter[0][fx] * temp_output[input_x][y];
            }
            output[x][y] = accum;
        }
    }
    return output;

}

// Convolve an image with a separable convolution kernel. Reused from my own code that I took last spring
//

int copy_edge(int val, int boundary)
{
    if (val<0)
        return 0;
    else if (val>=boundary)
        return boundary-1;
    else
        return val;

}
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
    SDoublePlane output(input.rows(), input.cols());


    for (int i=0; i< input.rows();i++) {
        for (int j = 0; j < input.cols(); j++) {
            double sum=0;
            for(int u=0;u<filter.rows();u++){
                for(int v=0;v<filter.cols();v++) {
                    int x=copy_edge(i+(u-filter.rows()/2),input.rows());
                    int y=copy_edge(j+(v-filter.cols()/2),input.cols());

                    sum+=filter[u][v]*input[x][y];
                }
            }
            output[i][j]=sum;

        }
    }
    // Convolution code here

    return output;
}

void template_convolution(const SDoublePlane &input,const SDoublePlane &template_image, vector<DetectedSymbol> &symbols,SDoublePlane &output, const Type t,int edge=0)
{

    for (int i=0; i< input.rows()-template_image.rows()+1;i++) {
        for (int j = 0; j < input.cols()-template_image.cols()+1; j++) {
            long double sum=0.0;
            for(int k=0;k<template_image.rows();k++){

                for(int l=0;l<template_image.cols();l++) {
                    int ik=floor((((input[i+k][j+l])/255.0))+0.5);
                    int tk=floor(((template_image[k][l])/255.0)+0.5);
                    sum+=ik*tk;
                    // cout<<sum;
                    if (edge==0)
                        sum+=((1-ik)*(1-tk));
                }
            }
            sum=(template_image.rows()*template_image.cols())-sum;
            output[i][j]=sum;
            if (sum<16)
            {
                DetectedSymbol s;
                s.row=i;
                s.col=j;
                s.width=template_image.cols();
                s.height=template_image.rows();

                s.type = t;
                s.confidence = 1.0;
                s.pitch = get_pitch(s.row, s.col, s.row+s.height-1, s.col+s.width-1);
                symbols.push_back(s);
            }

        }

    }

}

char get_pitch(int _top, int _left, int _bottom, int _right)
{
    int total_staffs=staffs.size()/5;

    for(int i=0;i<total_staffs;i++)
    {
        int beginning=staffs.at(i*5);
        int mid=(_top+_bottom)/2;
        if(_bottom<beginning-avg_space)
        {
            return 'B';

        }
        else if(_top>(beginning-avg_space) and _bottom<beginning)
            return 'G';
        else if(_top>(beginning) and _bottom<(beginning+avg_space))
            return 'E';
        else if (_top>(beginning+avg_space) and _bottom<(beginning+2*avg_space))
            return 'C';
        else if (_top>(beginning+2*avg_space) and _bottom<(beginning+3*avg_space))
            return 'A';
        else if (_top>(beginning+3*avg_space) and _bottom<(beginning+4*avg_space))
            return 'F';
        else if (_top>(beginning+4*avg_space) and _bottom<(beginning+5*avg_space))
            return 'D';
        else if (_top>(beginning+5*avg_space) and _bottom<(beginning+6*avg_space))
            return 'B';
        else
            return '-';
    }

    return '-';

}

void template_convolution_edgemap(const SDoublePlane &input,const SDoublePlane &template_image, vector<DetectedSymbol> &symbols,SDoublePlane &output, const Type t)
{

    for (int i=0; i< input.rows()-template_image.rows()+1;i++) {
        for (int j = 0; j < input.cols()-template_image.cols()+1; j++) {
            long double sum=0.0;
            for(int k=0;k<template_image.rows();k++){

                for(int l=0;l<template_image.cols();l++) {
                    sum+=input[i+k][j+l]*template_image[k][l];
                    // cout<<sum;

                }
            }
            output[i][j]=sum;
            if (sum<16)
            {
                DetectedSymbol s;
                s.row=i;
                s.col=j;
                s.width=template_image.cols();
                s.height=template_image.rows();
                //get_pitch(i,s.height);

                s.type = t;
                s.confidence = 1.0;
                s.pitch = get_pitch(s.row, s.col, s.row+s.height-1, s.col+s.width-1);
                symbols.push_back(s);
            }

        }

    }

}


