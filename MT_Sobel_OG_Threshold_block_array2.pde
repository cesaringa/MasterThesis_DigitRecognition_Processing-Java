int number_vectors = 8;              // Number of elements to represent the oriented gradients
int gridx = 2;                       // Number of divisions of oriented gradient matrix on X axis 
int gridy = 3;                       // Number of divisions of oriented gradient matrix on Y axis 
int NUMBER_SAMPLES = 20;              // Number of images that will be tested
int DIGITS = 2;          
PImage [][] matrix_digits = new PImage [DIGITS][NUMBER_SAMPLES];
PImage [][] grayscale = new PImage [DIGITS][NUMBER_SAMPLES];
PImage [][] threshold = new PImage [DIGITS][NUMBER_SAMPLES];

void setup() { 
  //size(1500, 1000);
  for (int j =0; j < DIGITS; j++) {
    for (int i=0; i < NUMBER_SAMPLES; i++) {
      matrix_digits[j][i] = loadImage(+j+"-"+i+".jpg");
      grayscale[j][i] = grayscale(matrix_digits[j][i] );
      threshold[j][i] = adaptive_threshold(grayscale[j][i] );    
    }
  }
    for (int k =0; k < DIGITS; k++) 
    for (int j =0; j < DIGITS; j++) {
    for (int i=0; i < NUMBER_SAMPLES; i++) {
      float [][] new_matrix_1 = feature_vector(threshold[k][0] );
      float [][] new_matrix_2 = feature_vector(threshold[j][i] );
      float c = 0.0;
      c = compare(new_matrix_1, new_matrix_2);
      
      println(""+k+"-0 and "+j+"-"+i+" =" ,c );
    }
  }
}

 PImage grayscale(PImage img) {
  PImage temp = createImage(img.width, img.height, RGB);
  temp.loadPixels();
  img.loadPixels();
  for (int y = 0; y < img.height; y++)
    for (int x = 0; x < img.width; x++) {
      int imgIndex = x + y * img.width;
      float r = red (img.pixels[imgIndex]);
      float g = green(img.pixels[imgIndex]);
      float b = blue(img.pixels[imgIndex]);      
      temp.pixels[imgIndex] = color(0.21*r + 0.72*g + 0.07*b); // load intensity
    }
  temp.updatePixels();
  return temp;
 }

//////////////////////////////////// FUNCTION FOR ADAPTIVE THRESHOLD ////////////////////////////////////
PImage adaptive_threshold(PImage img) {
  PImage temp = createImage(img.width, img.height, RGB);
  temp.loadPixels();
  img.loadPixels();
   
  float t = 0;                    // final threshold value
  float u1 = 0;                   // variable 1 for adaptive threshold
  float u2 = 0;                   // variable 2 for adaptive threshold
  //int n1 = 0;                     // Number of pixels group 1
  int n2 = 0;                     // Number of pixels group 2
  float t1=0;                     // Threshold value from 1st iteration
  float d = 0;                    // Difference between threshold values from iterations
  float sum = 0;                  // Sum up brightness from the whole image
  float sum1 = 0;                 // Sum up brightness from the first group of the image (<= threshold value)
  float sum2 = 0;                 // Sum up brightness from the second group of the image (> threshold value)
  
  for (int y = 0; y < img.height; y++)
      for (int x = 0; x < img.width; x++) {
        int imgIndex = x + y * img.width;
        sum += brightness(img.pixels[imgIndex]);
      }
    t = sum/(img.width * img.height);
  
  do {
      int n1 = 0;
      for (int y = 0; y < img.height; y++)
        for (int x = 0; x < img.width; x++) {
          int imgIndex = x + y * img.width;
          if (brightness(img.pixels[imgIndex]) <= t) {
            sum1 += brightness(img.pixels[imgIndex]);
            n1++;
          }
        }  
      n2 = ((img.width * img.height) - n1);
      sum2 = (sum - sum1);
      u1 = sum1 / n1;
      u2 = sum2 / n2;   
      t1 = (u1 + u2) / 2;
      d = abs(t1-t);
      t=t1;
    } while (d>5);
     
    for (int y = 0; y < img.height; y++)
    for (int x = 0; x < img.width; x++) {
      int imgIndex = x + y * img.width; 
     /* 
      float r = red (img.pixels[imgIndex]);
      float g = green (img.pixels[imgIndex]);
      float b = blue (img.pixels[imgIndex]);      
      float threshold_image = 0.21*r + 0.72*g + 0.07*b; //It is working when adaptive threshold is calculated with color intensity
*/
      float threshold_image = brightness(img.pixels[imgIndex]); // When extracting brightness from grayscale, threshold does not work properly
      
      if ( threshold_image <= t )
        threshold_image = 0;                  // black
      else threshold_image= 256;              // white
      
      temp.pixels[imgIndex] = color(threshold_image);
    }
  temp.updatePixels();
  return temp;
 }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
//////////////////////////////////// FUNCTION FOR GRADIENT IMAGE REPRENSATION ///////////////////////////////////////////////////////////////////
float [][] feature_vector(PImage img) { 
  
  float[][] sobelfilter_vert ={ { -1, 0, 1 }, //vertical sobel filter
                                { -2, 0, 2 }, 
                                { -1, 0, 1 } };
  float[][] sobelfilter_hori ={ { -1, -2, -1 }, //horizontal sobel filter
                                { 0, 0, 0 }, 
                                { 1, 2, 1 } };   
  float theta;
  float [][] new_matrix = new float [(img.width -2) * (img.height-2)][number_vectors];
  float [][] block_matrix = new float [gridx * gridy][number_vectors];
  
  ////////////////////////// Aplication of Sobel Filter ///////////////////////
  // Loop through every pixel in the image.
  for (int y = 1; y < img.height-1; y++) { // Skip top and bottom edges
    for (int x = 1; x < img.width-1; x++) { // Skip left and right edges
      float gx = 0, gy = 0; // Sobel sum for this pixel
      for (int ky = -1; ky <= 1; ky++) {
        for (int kx = -1; kx <= 1; kx++) {
          // Calculate the adjacent pixel for this Sobel point
          int index = (y + ky)*img.width + (x + kx);
    
          float pre_processed_img = brightness(img.pixels[index]);
   
          // Multiply adjacent pixels based on the Sobel filter values
          gx += sobelfilter_vert[ky+1][kx+1] * pre_processed_img;
          gy += sobelfilter_hori[ky+1][kx+1] * pre_processed_img;
        }
      }
      //////////////////////////////////////////////////////////////////////////////
      
      ////////////////////////// Gradient Image Representation ///////////////////////
      theta = atan2(gy, gx);
      if (theta < 0) theta += TWO_PI;
      
      // comparing the direction of the gradient in order to place into the matrix
      for (int i = 0; i < number_vectors; i++) {
        if ((theta>=TWO_PI/number_vectors*i) && (theta<TWO_PI/number_vectors*(i+1)))   
          new_matrix [(y-1)*(img.width-2) + (x-1)] [i] = sqrt(gx*gx + gy*gy);
        else 
          new_matrix [(y-1)*(img.width-2) + (x-1)] [i] = 0;
      }
      //////////////////////////////////////////////////////////////////////////////      
      // !!!!!
      // x y
      // p, q, ind
      int p = (int) (x * gridx/img.width);
      int q = (int) (y * gridy/img.height);
      
      int index_block = q*gridx+p;
      for (int j = 0; j < number_vectors; j++) {
      block_matrix[index_block][j] += new_matrix [(y-1)*(img.width-2) + (x-1)] [j];
      }
      
    }
  }
    float sum = 0;
    for (int i = 0; i < gridx * gridy; i++)
     for (int j = 0; j < number_vectors; j++) {
      sum += sq(block_matrix[i][j]);
      }

    for (int i = 0; i < gridx * gridy; i++)
     for (int j = 0; j < number_vectors; j++) {
      block_matrix[i][j] /= sum;
      }

  return block_matrix;
}

//////////////////////////////////// FUNCTION FOR COMPARING SIMILARITY OF 2 IMAGES ////////////////////////////////////
float compare(float [][] m1, float [][] m2) {
  float s12 = 0.0, s1 = 0.0, s2= 0.0, c=0; 
  for (int y = 0; y < gridx * gridy; y++) {
    for (int x = 0; x < number_vectors; x++) { 
    s12 += m1 [y][x] * m2 [y][x];
    s1 += sq(m1 [y][x]);
    s2 += sq(m2[y][x]);
    }
  }
  c = s12/(sqrt(s1)*sqrt(s2));
  return c;
} 
