/* ParticleExample
 no recenter on particle files
 */

#include <slsimlib.h>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>

#include "gridmap.h"
#include <math.h>
#include "utilities.h"
#include "utilities_slsim.h"
#include "lens_halos.h"
#include "disk.h"
#include "concave_hull.h"
#include "grid_maintenance.h"


using namespace std;

int main(int arg,char **argv){
  
  double range = 10*arcsecTOradians;
  double final_pixel_size = 0.1*arcsecTOradians;
  int over_sample_factor = 4;  // oversapling of PSF
  
  double mag = 23;  // magnitude of source
  double Reff = 0.2;  // R effective of source in arcseconds
  double pos_angle = PI/4;  // position angle in radians
  double nn = 4; // Sersic index
  double fratio = 0.5; // axis ratio
  
  //double count_zero_point = 25.7;  // magnitude that produces 1 e- per second
  double count_zero_point = 15.7;  // magnitude that produces 1 e- per second
  double background_rms_e = 0.003;  // rms background noise in e- per second per pixel
  const std::vector<double> exp_times = {560.52,560.52,560.52,560.52,89.52,89.52};  // exposure times of dithers in seconds
  Point_2d center(0,0);  // center of image
 
  Utilities::RandomNumbers_NR noise_ran;  // random number generator
  
  std::map<Band,double> zeropoints;
  zeropoints.emplace(Band::EUC_VIS,25.); // converted to AB mag -> e-/s
  zeropoints.emplace(Band::EUC_Y,29.8);// these are from AB mag -> e-/s
  zeropoints.emplace(Band::EUC_J,30.0);
  zeropoints.emplace(Band::EUC_H,30.0);
  
  // you can set the zeropints like this all at once or when you make each source
  SourceColored::setMagZeroPoints(zeropoints);
  
// **************************************************************************
  
  int Npixels = (int)(range/final_pixel_size);  // number of pixels for the final image
  
  SourceSersic source(mag,Reff,pos_angle, nn, fratio, 1,Band::EUC_VIS);
  
  // this is how you might set and change the band
  //source.setMag(mag,Band::EUC_VIS, zeropoints[Band::EUC_Y]);
  //source.setActiveBand(Band::EUC_Y);
 
  ObsVIS obs(Npixels,Npixels
              ,exp_times
              ,over_sample_factor
              ,final_pixel_size
              ,background_rms_e
              ,Utilities::vec_sum(exp_times)
              );
  
  
  // By default the psf is taken to be Gaussian, but a psf map can be input
  //obs.setPSF("psf_file.fits");
  
  // if no PSF map is provided Gaussian smoothing can be done by setting the FWHM in arcseconds
  obs.setSeeing(1.5*final_pixel_size/arcsecTOradians);
  obs.setZeropoint( zeropoints[Band::EUC_VIS] );
  
  PixelMap<float> oversampled_image(center.x,over_sample_factor*Npixels,final_pixel_size/over_sample_factor);
  // or PixelMap<float> oversamples_image(center.x,obs.getNxInput(),obs.getPixelSize()/over_sample_factor
  //                               ,PixelMapUnits::count_per_sec);
  
  oversampled_image.AddSource(source);  // No lensing here, just add the source.
  
  PixelMap<float> output_image(center.x,Npixels,final_pixel_size
                                 ,PixelMapUnits::count_per_sec);
  
  PixelMap<float> error_maps(output_image);  // error map with same dimensions
 
  // this will apply the psf, add noise and resample the image
  obs.Convert<float>(oversampled_image
                      ,output_image
                      ,error_maps
                      ,true     // apply psf
                      ,true     // include noise
                      ,noise_ran
                      );
  
  
  oversampled_image.printFITS("imput_image.fits");
  output_image.printFITS("output_image.fits");
}
