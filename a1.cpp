#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <DrawText.h>
#include <sstream>
#include <time.h>


#include "convolution.cpp"


#define TWO_PI 6.2831853071795864769252866
#define canny_threshold 1


using namespace std;


// DetectedSymbol class may be helpful!
//  Feel free to modify.
//

class Point
{
public:
    int x,y;
    Point(){x=0;y=0;}
    Point(int u,int v){x=u;y=v;}

};
class Stave
{
public:
    vector<int> height;

};



void write_image(const string &filename, const SDoublePlane &output, int grayscale,int normalize=0);
SDoublePlane canny_edge_detector(const SDoublePlane &input,int region_size=5);
SDoublePlane apply_gaussian(const SDoublePlane &input, double sigma);
SDoublePlane image_rescale(const SDoublePlane &input_image);
//void template_convolution(const SDoublePlane &input,const SDoublePlane &template_image, vector<DetectedSymbol> &symbols,SDoublePlane &output, const Type t,int edge)
//;
//void template_convolution_edgemap(const SDoublePlane &input,const SDoublePlane &template_image, vector<DetectedSymbol> &symbols,SDoublePlane &output, const Type t);


// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc.
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose.

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
	int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

	// if any of the coordinates are out-of-bounds, truncate them
	top = min( max( top, 0 ), input.rows()-1);
	bottom = min( max( bottom, 0 ), input.rows()-1);
	left = min( max( left, 0 ), input.cols()-1);
	right = min( max( right, 0 ), input.cols()-1);

	// draw top and bottom lines
	for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
	// draw left and right lines
	for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;

  }
}


// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<struct DetectedSymbol> &symbols)
{
  ofstream ofs(filename.c_str());

  for(int i=0; i<symbols.size(); i++)
	{
	  const DetectedSymbol &s = symbols[i];
	  ofs << s.row << " " << s.col << " " << s.width << " " << s.height << " ";
	  if(s.type == NOTEHEAD)
	ofs << "filled_note " << s.pitch;
	  else if(s.type == EIGHTHREST)
	ofs << "eighth_rest _";
	  else
	ofs << "quarter_rest _";
	  ofs << " " << s.confidence << endl;
	}

}

// Function that outputs a visualization of detected symbols
void  write_detection_image(const string &filename, const vector<DetectedSymbol> &symbols, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];
  for(int i=0; i<3; i++)
	output_planes[i] = input;

  for(int i=0; i<symbols.size(); i++)
	{
	  const DetectedSymbol &s = symbols[i];

	  overlay_rectangle(output_planes[s.type], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 255, 2);
	  overlay_rectangle(output_planes[(s.type+1) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);
	  overlay_rectangle(output_planes[(s.type+2) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);

	  if(s.type == NOTEHEAD)
	{
	  char str[] = {s.pitch, 0};
	  draw_text(output_planes[0], str, s.row, s.col+s.width+1, 0, 2);
	  draw_text(output_planes[1], str, s.row, s.col+s.width+1, 0, 2);
	  draw_text(output_planes[2], str, s.row, s.col+s.width+1, 0, 2);
	}
	}

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}



// The rest of these functions are incomplete. These are just suggestions to
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a separable convolution kernel
//

// Apply a sobel operator to an image, returns the result
//
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
    SDoublePlane output(input.rows(), input.cols());
    SDoublePlane row_filter(1,3);
    SDoublePlane col_filter(3,1);



    //Sobel Filter
    for(int i=0;i<3;i++)
    {
        if(_gx) {
            //Horizontal Edges
            row_filter[0][i] = i < 2 ? (i + 1 )/8.0: 1/8.0;
            col_filter[i][0] = ( i-1)/8.0;

        }
        else
        {
            // Vertical Edges
            col_filter[i][0] = i < 2 ? (i + 1)/8.0 : 1/8.0;
            row_filter[0][i] = (i-1)/8.0;

        }
    }


    // Implement a sobel gradient estimation filter with 1-d filters

    output=  convolve_separable(input, row_filter,col_filter);
    return output;
}

// Apply an edge detector to an image, returns the binary edge map
//
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
    SDoublePlane output(input.rows(), input.cols());
    SDoublePlane vertical_output(input.rows(), input.cols());
    SDoublePlane horizontal_output(input.rows(), input.cols());
    SDoublePlane edge_output(input.rows(), input.cols());

  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold
    vertical_output=sobel_gradient_filter(input,1);
    horizontal_output=sobel_gradient_filter(input,0);

    write_image("test_images/horizontal_edgs.png",horizontal_output,1);
    write_image("test_images/vertical_edges.png",vertical_output,1);

    for(int x=0;x<input.rows();x++)
    {
        for(int y=0;y<input.cols();y++)
        {
            double edge_strength=sqrt((vertical_output[x][y]*vertical_output[x][y])+(horizontal_output[x][y]*horizontal_output[x][y]));

            //use threshold to produce binary image
            if(edge_strength > thresh)
                output[x][y]=1.0;
            else
                output[x][y]=0;

        }
    }

    write_image("sobel.png",output,1);

  return output;
}




void calculate_distance_map_chamfer(const SDoublePlane &edge_input, SDoublePlane &distance_map)
{

    for(int i=0;i<edge_input.rows();i++)
    {
        for(int j=0;j<edge_input.cols();j++)
        {

            if (edge_input[i][j]==1)
            {
                distance_map[i][j]=0;
              //  cout<<distance_map[i][j]<<endl;
            }
            else
            {
                distance_map[i][j]=1000;
               // cout<<distance_map[i][j]<<endl;
            }


        }



    }



    double d1=0,d2=0,d3=0,d4=0;
   // cout<<edge_input.cols()<<endl;
    for(int j=0;j<edge_input.cols();j++)
    {
        for(int i=0;i<edge_input.rows();i++)
        {
         //   cout<<i<<j<<endl;
            if (distance_map[i][j]>0)
            {
                d1=i-1>0?1+distance_map[i-1][j]:distance_map[i][j];
                d2=i-1>0 && j-1>0? sqrt(2)+distance_map[i-1][j-1]:distance_map[i][j];
                d3=j-1>0?1+distance_map[i][j-1]:distance_map[i][j];
                d4=i+1<edge_input.rows() && j-1>0?sqrt(2)+distance_map[i+1][j-1]:distance_map[i][j];

                distance_map[i][j]=std::min(d1, std::min(d2, std::min(d3, d4)));
                //cout<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<endl;
              //  cout<<distance_map[i][j]<<endl;
            }

        }

    }


    for(int j=edge_input.cols()-1;j>=0;j--)
    {
        for(int i=edge_input.rows()-1;i>=0;i--)
        {

            if (distance_map[i][j]>0)
            {
                d1=i+1<edge_input.rows()?1+distance_map[i+1][j]:distance_map[i][j];
                d2=i+1<edge_input.rows() && j+1<edge_input.cols()? sqrt(2)+distance_map[i+1][j+1]:distance_map[i][j];
                d3=j+1<edge_input.cols()?1+distance_map[i][j+1]:distance_map[i][j];
                d4=i-1>0 && j+1<edge_input.cols()?sqrt(2)+distance_map[i-1][j+1]:distance_map[i][j];

                distance_map[i][j]=std::min(distance_map[i][j],std::min(d1, std::min(d2, std::min(d3, d4))));
            }

        }

    }





}



void template_matching_edge_map(const SDoublePlane &input)
{
    double sigma=1.4;
   SDoublePlane edge_input=find_edges(input,8);//apply_gaussian(input,sigma);
 //  SDoublePlane edge_input=apply_gaussian(input,sigma);
   // edge_input=find_edges(edge_input,7.0);

    write_image("edges.png",edge_input,1,1);
 //   edge_input=canny_edge_detector(edge_input,5);
   // write_image("canny_edges.png",edge_input,1,1);


    SDoublePlane output(input.rows(), input.cols());

    SDoublePlane distance_map(edge_input.rows(),edge_input.cols());
   // calculate_distance_map(edge_input,distance_map);
    calculate_distance_map_chamfer(edge_input,distance_map);

    cout<<"Here"<<endl;
    vector<DetectedSymbol> symbols;
    SDoublePlane template_image;

    for(int i=0;i<3;i++){
        std::string s;
        std::stringstream out;
        out << (i+1);
        s = out.str();
        string filename="template"+s+".png";
        //cout<<filename<<endl;
        template_image= SImageIO::read_png_file(filename.c_str());
     //   template_image=image_rescale(template_image);
        SDoublePlane template_input=find_edges(template_image,6.0+i*i);//apply_gaussian(template_image,sigma);
       // SDoublePlane template_input=apply_gaussian(template_image,sigma);
     //   template_input=find_edges(template_input,7.0);
       // template_input=canny_edge_detector(template_input,2);

        write_image("template_canny_"+s+".png",template_input,1,1);

        template_convolution_edgemap(distance_map,template_input,symbols,output,(Type)i);
        cout<<"Total Symbols Detected: "<<symbols.size()<<endl;
    }


    // Convolution code here
    write_detection_txt("detected5.txt", symbols);
    write_detection_image("detected5.png", symbols, input);
    write_image("scores5.png",output,1);




}


void write_image(const string &filename, const SDoublePlane &output, int grayscale=0,int normalize)
{
    SDoublePlane normalized(output.rows(),output.cols());
    if(normalize==1)
    {
        for(int i=0;i<output.rows();i++)
        {
            for (int j=0;j<output.cols();j++)
            {
                if(output[i][j]==1)
                    normalized[i][j]=255;
            }
        }
        SImageIO::write_png_file(filename.c_str(), normalized,normalized,normalized);
    }

    else
    {
        if (grayscale==1)
            SImageIO::write_png_file(filename.c_str(), output,output,output);

    }


}
double calculate_Gaussian(int x, int y, int sigma)
{
    return ((1/(TWO_PI*(sigma*sigma)))*exp(-(((double)(x*x)+(y*y))/(2*(sigma*sigma)))));

}

SDoublePlane apply_gaussian(const SDoublePlane &input, double sigma)
{



    int filterSize = 5;
    SDoublePlane filter(filterSize, filterSize) ;
    double value=0,value2=0.0;
    for(int y = 0; y < filterSize; y++) {
        for(int x = 0; x < filterSize; x++)
        {
            filter[x][y] = calculate_Gaussian(x-(filterSize/2),y-(filterSize/2),sigma);
            value+=filter[x][y];

        }
        // cout<<endl;
    }

    for(int y = 0; y < filterSize; y++) {
        for(int x = 0; x < filterSize; x++)
        {
           filter[x][y] = filter[x][y]/value;
            // value2+=H(x,y,0,0);

        }
        // cout<<endl;
    }



    SDoublePlane output_image = convolve_general(input, filter);
    write_image("gaussian.png",output_image,1);
    return output_image;


}


void apply_separable(const SDoublePlane &input, int filter_size)
{


    SDoublePlane row_filter(1,filter_size);
    SDoublePlane col_filter(filter_size,1);

    //Mean Filter
    for(int i=0;i<filter_size;i++)
    {
        row_filter[0][i]=1/(double)filter_size;
        col_filter[i][0]=1/(double)filter_size;
    }

    cout<<"Seperable Filter Created with filter size "<<row_filter.rows()<<", Output File: seperable.png "<<endl;
    SDoublePlane output_image = convolve_separable(input, row_filter,col_filter);
    write_image("seperable.png",output_image,1);


}





void template_matching(const SDoublePlane &input,int p7=0 )
{

    SDoublePlane output(input.rows(), input.cols());
    vector<DetectedSymbol> symbols;
    SDoublePlane template_image;

    for(int i=0;i<3;i++){
        std::string s;
        std::stringstream out;
        out << (i+1);
        s = out.str();
        string filename="template"+s+".png";
        //cout<<filename<<endl;
        template_image=SImageIO::read_png_file(filename.c_str());
      //  template_image=image_rescale(template_image);
        template_convolution(input,template_image,symbols,output,(Type)i);
        cout<<"Total Symbols Detected (General Matching): "<<symbols.size()<<endl;
    }


    // Convolution code here
    write_detection_txt("detected.txt", symbols);
    write_detection_image("detected4.png", symbols, input);
    write_image("scores4.png",output,1);

    if(p7==1)
    {
        write_detection_txt("detected7.txt", symbols);
        write_detection_image("detected7.png", symbols, input);
    }

   // return output;
}

//
// This main file just outputs a few test images. You'll want to change it to do
//  something more interesting!
//

SDoublePlane canny_edge_detector(const SDoublePlane &input, int region_size)
{

  //  int region_size=5;
    SDoublePlane working_input=input;
    SDoublePlane output(input.rows(),input.cols());

    int max=1;
    for(int x=0;x<working_input.rows();x++)
    {
        for(int y=0;y<working_input.cols();y++)
        {
            max=1;
            for(int a=0;a<region_size;a++)
            {
                for(int b=0;b<region_size;b++)
                {
                    int i=copy_edge(x+(a-region_size/2),input.rows());
                    int j=copy_edge(y+(b-region_size/2),input.cols());
                    if(working_input[i][j]>working_input[x][y]) {
                        max = 0;
                        break;
                    }
                }
            }
            if(max==0)
            {
                working_input[x][y]=0;

            }

            if(working_input[x][y]>canny_threshold)
            {
                working_input[x][y]=1;
                output[x][y]=255;

            }
            else
            {
                working_input[x][y]=0;
                output[x][y]=0;
            }
            //use threshold to produce binary image


        }
    }
    write_image("canny_output.png",output,1);
    return working_input;

}

void hough_lines(const SDoublePlane &input_image)
{

    SDoublePlane input(input_image.rows(), input_image.cols());


    double sigma=1.4;
    SDoublePlane working_input=apply_gaussian(input_image,sigma);
    working_input=sobel_gradient_filter(input_image,0);
    input=canny_edge_detector(working_input);//cann=sobel_gradient_filter(input_image,0);


    int accumulator[input.rows()];
    memset( accumulator, 0, input.rows()*sizeof(int) );


    for(int u=0;u<input.rows();u++)
    {
        for(int v=0;v<input.cols();v++)
        {
            if(input[u][v]==1)
            {
                //double x=(double)u-0;
                accumulator[u]+=1;

            }
        }

    }

    SDoublePlane output_planes[2];
    SDoublePlane blue_plane=input_image;
    for(int i=0; i<2; i++)
        output_planes[i] = input_image;

    double avg_height=0;
    int prev_u=0,five=0,staff_index=0;
    int total_staffs=1;
    for(int u=0;u<input.rows();u++)
    {
        if(accumulator[u]>input.cols()/5)
        {

            //cout<<u<<accumulator[u]<<endl;
            if(u-prev_u>2 && five<5)
            {
                //cout<<u<<" "<<(u-prev_u)<<endl;
                staffs.push_back(u);
                avg_height+=(five==0?0:(u-prev_u));
               // cout<<avg_height<<endl;
                prev_u=u;
                five++;
                staff_index++;
            }

            for(int v=0;v<input.cols();v++) {
                blue_plane[u][v]=255;
                output_planes[0][u][v]=0;
                output_planes[1][u][v]=0;
            }


            if(five==5)
                total_staffs++;
            five=five==5?0:five;

        }


    }
    if(avg_height>0)
        avg_space=floor(avg_height/((total_staffs-1)*4));
    //cout<<avg_height/8<<endl;
    cout<<"Total Lines Detected: "<<(total_staffs-1)*5<<endl;
    SImageIO::write_png_file("staves.png",output_planes[0],output_planes[1],blue_plane);



}


SDoublePlane image_rescale(const SDoublePlane &input_image)
{
    double rescaling_factor=(double)avg_space/11; //11 for template 1
    int modified_rows=ceil(input_image.rows()/rescaling_factor);
    int modified_cols=ceil(input_image.cols()/rescaling_factor);
    SDoublePlane rescaled_image(modified_rows,modified_cols);

    for(int x=0;x<modified_rows;x++) {
        for (int y = 0; y < modified_cols; y++) {
            int new_x=floor(x*rescaling_factor);
            int new_y=floor(y*rescaling_factor);
            if(new_x<input_image.rows() && new_y<input_image.cols())
                rescaled_image[x][y]=input_image[new_x][new_y];
            else
                rescaled_image[x][y]=128;
        }
    }


    write_image("rescaled.png",rescaled_image,1);
    return rescaled_image;


    /*double rescaling_factor=(double)avg_space/11; //11 for template 1
    int modified_rows=ceil(input_image.rows()/rescaling_factor);
    int modified_cols=ceil(input_image.cols()/rescaling_factor);
    int row_delete=modified_rows>=input_image.rows()?modified_rows-input_image.rows():input_image.rows()-modified_rows;
    int col_delete=modified_cols>=input_image.cols()?modified_cols-input_image.cols():input_image.cols()-modified_cols;
    cout<<modified_rows<<endl;
    cout<<modified_cols<<endl;
    SDoublePlane rescaled_image(modified_rows,modified_cols);
    int u=0,v=0;
    double p=0;
    double probability=abs(1-rescaling_factor);
    srand (time(NULL));

    for(int x=0;x<input_image.rows();x++)
    {
        for(int y=0;y<input_image.cols();y++)
        {
            int random=rand()%2;

            if(random==0 && row_delete>0 && p<probability )
            {
                rescaled_image[u][v]=

            }
        }





    }
    return rescaled_image;*/
}

int main(int argc, char *argv[])
{
    if(!(argc >= 2))
    {
        cerr << "usage: " << argv[0] << " input_image"<< "part no (Optional, Default:p7) " << endl;
        return 1;
    }
    string input_filename(argv[1]);
    string part_no("");
    if(argc>=3)
    {
        part_no=argv[2];
    }


    SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
    // test step 2 by applying mean filters to the input image

    if(part_no=="p2")
    {
        int filter_size=5;

        SDoublePlane mean_filter(filter_size,filter_size);
        for(int i=0; i<filter_size; i++)
            for(int j=0; j<filter_size; j++)
                mean_filter[i][j] = 1/(double)(filter_size*filter_size);

        //part 2

        SDoublePlane output_image = convolve_general(input_image, mean_filter);
        write_image("general_convolution.png",output_image,1);
        cout<<" General Convolution Done. Output FIle: general_convolution.png, Filter Size 5"<<endl;
    }




    //part 3
    if(part_no=="p3")
        apply_separable(input_image,5);


    // Part 4
    if(part_no=="p4")
    {
        cout<<"Part 4: Template Matching "<<endl;
        template_matching(input_image);
    }


    //Part 5
   // find_edges(input_image,8.0);
    if(part_no=="p5")
    {
        cout<<"Part 5: Template Matching alternative approach "<<endl;
        template_matching_edge_map(input_image);

    }


    //Part6
   // SDoublePlane sobel=find_edges(input_image,8.0);
    if(part_no=="p6") {
        cout<<"Part 6: Hough Lines"<<endl;
        hough_lines(input_image);
    }

    if(part_no.length()==0)
    {
        //Part7
        cout<<"Part 7 Started"<<endl;
        hough_lines(input_image);
        input_image=image_rescale(input_image);
        template_matching(input_image,1);
        template_matching_edge_map(input_image);

    }



   // calculate_distance_map_chamfer(mean_filter,distance);
    //output_image= template_ma;//find_edges(input_image,8.0);
    //write_image("sobel.png",output_image,1);


//   apply_separable(input_image,9);




      // randomly generate some detected symbols -- you'll want to replace this
      //  with your symbol detection code obviously!
//


  /*  for(int i=0; i<10; i++)
    {
      DetectedSymbol s;
      s.row = rand() % input_image.rows();
      s.col = rand() % input_image.cols();
      s.width = 20;
      s.height = 20;
      s.type = (Type) (rand() % 3);
      s.confidence = rand();
      s.pitch = (rand() % 7) + 'A';
      symbols.push_back(s);
    }*/

//    write_detection_txt("detected.txt", symbols);
 //   write_detection_image("detected.png", symbols, input_image);
}
