#pragma once

#include <algorithm>    // std::transform
#include <array>
#include <bitset>
#include <cmath>
#include <cstdarg>		// va_list, va_arg, va_start, va_end, va_copy
#include <cstdio>		// printf, sprintf,  
#include <cstdlib>		// rand, srand
#include <ctime>		// clock(), time()
#include <fcntl.h>
#include <fstream>
#include <functional>	// std::multiplies, std::plus, std::function, std::negate
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <new>			 
#include <numeric>		// inner_product, partial_sum, adjacent_difference, accumulate
//#include <omp.h>		// OpenMP
#include <process.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include "sys/types.h"	// stat f
#include "sys/stat.h"	// stat functions
#include <tuple>
#include <typeinfo>		//operator typeid
#include <type_traits>	// is_pod
//#include <unistd.h>
#include <utility>		// for std::move
#include <vector>
#if defined(_WIN32) || defined(_WIN64)
	#include <windows.h>
	#include "Shlwapi.h"	// Requires: Shlwapi.lib, Shlwapi.h, Shlwapi.dll (version 4.71 or later)
#endif
/***************************************************************************************************************************************************************************/
/************************************************************** Global typedef and namespace usage definitions *************************************************************/
/***************************************************************************************************************************************************************************/
//using namespace std::placeholders; 
//using namespace std;
using std::cout;
using std::endl;
typedef unsigned long long ULL;
typedef unsigned int uint;
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------- WED calculations parameters ----------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
enum DISK_WRITE_MODE { TEXT, BINARY };										// Experimental or simulated data
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------- WED calculations parameters ----------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
#define ANGLE_TO_RADIANS		( PI/180.0 )							// Convertion from angle to radians
#define PI						( 4 * atan( 1.0 ) )						// 4*pi/4 radians =   pi   radians = 180 degrees
#define START					true									// Used as an alias for true for starting timer
#define STOP					false									// Used as an alias for false for stopping timer
#define RIGHT					1										// Specifies that moving right corresponds with an increase in x position, used in voxel walk 
#define LEFT					-1										// Specifies that moving left corresponds with a decrease in x position, used in voxel walk 
#define UP						1										// Specifies that moving up corresponds with an increase in y/z position, used in voxel walk 
#define DOWN					-1										// Specifies that moving down corresponds with a decrease in y/z position, used in voxel walk 
#define X_INCREASING_DIRECTION	RIGHT									// [#] specifies direction (LEFT/RIGHT) along x-axis in which voxel #s increase
#define Y_INCREASING_DIRECTION	DOWN									// [#] specifies direction (UP/DOWN) along y-axis in which voxel #s increase
#define Z_INCREASING_DIRECTION	DOWN									// [#] specifies direction (UP/DOWN) along z-axis in which voxel #s increase
#define CONSOLE_WINDOW_WIDTH	80
#define STDOUT_2_DISK			false
#define OVERWRITE_RESULTS_OK	true
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------- WED calculations parameters ----------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
#define WED_TARGET_IMAGE_WIDTH ( 192 ) //mm
#define WED_TARGET_IMAGE_HEIGHT ( 192 ) //mm
#define WED_TARGET_IMAGE_THICKNESS ( 160 ) //mm
#define WED_TARGET_COLUMNS 1024
#define WED_TARGET_ROWS 1024
#define WED_TARGET_SLICES 128
#define WED_TARGET_VOXELS ( WED_TARGET_COLUMNS * WED_TARGET_ROWS * WED_TARGET_SLICES )
#define WED_TARGET_VOXEL_WIDTH ( 0.1875 ) // mm
#define WED_TARGET_VOXEL_HEIGHT ( 0.1875 ) // mm
#define WED_TARGET_VOXEL_THICKNESS ( 1.25) // mm
#define WED_TARGET_THRESHOLD_RSP 0.01 
#define WED_TARGET_X_ZERO_VOXEL_CENTER -90		// [mm] Coordinate in x direction corresponding to voxel 0 -96 -96.1875
#define WED_TARGET_Y_ZERO_VOXEL_CENTER 90			// [mm] Coordinate in y direction corresponding to voxel 0 96 96.1875
#define WED_TARGET_Z_ZERO_VOXEL_CENTER -6.875		// [mm] Coordinate in z direction corresponding to voxel 0 -6.875 -5.6250 -166.875 -168.1250
//#define WED_TARGET_X_ZERO_COORDINATE (WED_TARGET_X_ZERO_VOXEL_CENTER - WED_TARGET_VOXEL_WIDTH / 2)		// [mm] Coordinate in x direction corresponding to voxel 0 -96 -96.1875
//#define WED_TARGET_Y_ZERO_COORDINATE (WED_TARGET_Y_ZERO_VOXEL_CENTER - WED_TARGET_VOXEL_HEIGHT / 2)			// [mm] Coordinate in y direction corresponding to voxel 0 96 96.1875
//#define WED_TARGET_Z_ZERO_COORDINATE (WED_TARGET_Z_ZERO_VOXEL_CENTER - WED_TARGET_VOXEL_THICKNESS / 2)		// [mm] Coordinate in z direction corresponding to voxel 0 -6.875 -5.6250 -166.875 -168.1250
//#define WED_TARGET_X_ZERO_COORDINATE (WED_TARGET_X_ZERO_VOXEL_CENTER - WED_TARGET_VOXEL_WIDTH / 2)		// [mm] Coordinate in x direction corresponding to voxel 0 -96 -96.1875
//#define WED_TARGET_Y_ZERO_COORDINATE (WED_TARGET_Y_ZERO_VOXEL_CENTER - WED_TARGET_VOXEL_HEIGHT / 2)			// [mm] Coordinate in y direction corresponding to voxel 0 96 96.1875
//#define WED_TARGET_Z_ZERO_COORDINATE (WED_TARGET_Z_ZERO_VOXEL_CENTER - WED_TARGET_VOXEL_THICKNESS / 2)		// [mm] Coordinate in z direction corresponding to voxel 0 -6.875 -5.6250 -166.875 -168.1250
#define WED_TARGET_X_ZERO_COORDINATE -90		// [mm] Coordinate in x direction corresponding to voxel 0 -96 -96.1875
#define WED_TARGET_Y_ZERO_COORDINATE 90			// [mm] Coordinate in y direction corresponding to voxel 0 96 96.1875
#define WED_TARGET_Z_ZERO_COORDINATE -6.875		// [mm] Coordinate in z direction corresponding to voxel 0 -6.875 -5.6250 -166.875 -168.1250
//#define WED_TARGET_X_ZERO_COORDINATE -96		// [mm] Coordinate in x direction corresponding to voxel 0 -96 -96.1875
//#define WED_TARGET_Y_ZERO_COORDINATE 96			// [mm] Coordinate in y direction corresponding to voxel 0 96 96.1875
//#define WED_TARGET_Z_ZERO_COORDINATE -6.875		// [mm] Coordinate in z direction corresponding to voxel 0 -6.875 -5.6250 -166.875 -168.1250
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------- WED calculations variables -----------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//"C:\\Users\\Blake\\Documents\\Visual Studio 2010\\Projects\\robust_pct\\robust_pct"
char* current_directory, *targets_directory, *WED_results_directory, *target_volume_directory;
char targets_folder[] = "bap_coordinates";
char target_volume_folder[] = "RStP_DICOM_PHANTOM";
char WED_results_folder[] = "WED_Results";
char WED_results_basename[] = "WED_Results";
char target_volume_filename[] = "RStP_Phantom";
const char targets_base_name[] = "proxi_distal";	

char EXECUTION_DATE[9];
unsigned int NUM_RUN_ARGUMENTS;
char** RUN_ARGUMENTS;
bool DATA_PATH_PASSED;
bool DEBUG_ON = false;

float* RSP_Phantom_image_h, *RSP_Phantom_image_d;
int num_targets;
double* target_x_h, * target_x_d;
double* target_y_h, * target_y_d;
double* target_z_h, * target_z_d;
double* WED_results;
std::vector<int> entry_voxel_x;
std::vector<int> entry_voxel_y;
std::vector<int> entry_voxel_z;
std::vector<int> beam_angles;
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------- Execution timer variables ------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
clock_t program_start, program_end, pause_cycles = 0;