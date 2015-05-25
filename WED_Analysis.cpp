#pragma once
/***********************************************************************************************************************************************************************************************************************/
/********************************************************************************** Proton CT Preprocessing and Image Reconstruction Code ******************************************************************************/
/***********************************************************************************************************************************************************************************************************************/
#include "WED_Analysis.h"

void apply_execution_arguments(unsigned int, char**);

int calculate_voxel( double, double, double );
double distance_remaining( double, double, int, int, double, int );
double edge_coordinate( double, int, double, int, int );
double path_projection( double, double, double, int, double, int, int );
double corresponding_coordinate( double m, double x, double x0, double y0 );

void get_dir_filenames(std::vector<std::string> &, const std::string &);
void get_dir_filenames_matching(const std::string &, const std::string &, std::vector<std::string> &, std::vector<std::string> &);
std::string terminal_response(char*) ;
void timer( bool, clock_t, clock_t);
void exit_program_if(bool);

void print_section_separator(char );
void construct_header_line( char*, char, char* );
void print_section_header( char*, char );
void print_section_exit( char*, char* );
char((&current_MMDD( char(&)[5]))[5]);
char((&current_MMDDYYYY( char(&)[9]))[9]);
void set_execution_date();
void set_IO_paths();

void WED_analysis();
void read_Phantom();
void read_txt_Phantom();
void read_bin_Phantom();
void read_slice( int);
void write_Phantom();
void write_slice( char*, const char*, double* &, int, int, int );
void combine_image_slices(const char*, char*, char*, uint, uint, uint, uint, DISK_WRITE_MODE );
void find_target_beam_angles();
void read_target_coordinates( int );
void print_targets();
double*& return_WED_results_old( int);
double*& return_WED_results( int);
void WED_take_2D_step
( 
	const int, const int, const int, const double, const double, const double, const double, const double, const double, 
	const double, const double, const double, double&, double&, double&, int&, int&, int&, int&, double&, double&, double&, double&	
);
void WED_take_3D_step
( 
	const int, const int, const int, const double, const double, const double, const double, const double, const double, 
	const double, const double, const double, double&, double&, double&, int&, int&, int&, int&, double&, double&, double&, double&	
);
double calculate_target_WED_old( float, float, float, float );
double calculate_target_WED( double, double, double, double, double );
void write_phantom_entries(int);
void write_WED_results( double*&, int);
/***********************************************************************************************************************************************************************************************************************/
/********************************************************************************** Proton CT Preprocessing and Image Reconstruction Code ******************************************************************************/
/***********************************************************************************************************************************************************************************************************************/
int main(unsigned int num_arguments, char** arguments)
{
	timer( START, program_start, program_end );
	

	set_IO_paths();
	WED_analysis();
	printf("%s\n", current_directory );
	exit_program_if(true);
}
void apply_execution_arguments(unsigned int num_arguments, char** arguments)
{
	NUM_RUN_ARGUMENTS = num_arguments;
	RUN_ARGUMENTS = arguments; 

	if( NUM_RUN_ARGUMENTS > 1)
		DATA_PATH_PASSED = true;
	else
		DATA_PATH_PASSED = false;

	
	// n =				  1			   1		   2			2		   3		   3	 ...	 n			 n       
	// i =	 0			  1			   2		   3			4		   5		   6	 ...   2n - 1  		 2n		  2n + 1
	// [program name][parameter 1][new val 1][parameter 2][new val 2][parameter 3][new val 3]...[parameter n][new val n][cfg path]
	//"C:\Users\Blake\Documents\Visual Studio 2010\Projects\robust_pct\robust_pct\settings.cfg"
	//"C:\Users\Blake\Documents\pCT_Data\object_name\Experimental\MMDDYYYY\run_number\Output\MMDDYYYY\Reconstruction\MMDDYYYY\settings.cfg"
	if( DATA_PATH_PASSED )
	{
		current_directory = (char*) calloc( strlen(RUN_ARGUMENTS[1])+1, sizeof(char));
		std::copy( RUN_ARGUMENTS[1], &RUN_ARGUMENTS[1][strlen(RUN_ARGUMENTS[1])], current_directory );
		print_section_header("Config file location passed as command line argument and set to : ",'*');
		print_section_separator('-');
		printf("%s\n", current_directory );
		print_section_separator('-');
	}
}
/***********************************************************************************************************************************************************************************************************************/
/********************************************************************************* Image Position/Voxel Calculation Functions (Host) ***********************************************************************************/
/***********************************************************************************************************************************************************************************************************************/
int calculate_voxel( double zero_coordinate, double current_position, double voxel_size )
{
	return abs( current_position - zero_coordinate ) / voxel_size;
}
double distance_remaining( double zero_coordinate, double current_position, int increasing_direction, int step_direction, double voxel_size, int current_voxel )
{
	/* Determine distance from current position to the next voxel edge.  path_projection is used to determine next intersected voxel, but it is possible for two edges to have the same distance in 
	// a particular direction if the path passes through a corner of a voxel.  In this case, we need to advance voxels in two directions simultaneously and to avoid if/else branches
	// to handle every possibility, we simply advance one of the voxel numbers and pass the assumed current_voxel to this function.  Under normal circumstances, this function simply return the 
	// distance to the next edge in a particual direction.  If the path passed through a corner, then this function will return 0 so we will know the voxel needs to be advanced in this direction too.
	*/
	int next_voxel = current_voxel + increasing_direction * step_direction;//  vz = 0, i = -1, s = 1 	
	double next_edge = edge_coordinate( zero_coordinate, next_voxel, voxel_size, increasing_direction, step_direction );
	return abs( next_edge - current_position );
}
double edge_coordinate( double zero_coordinate, int voxel_entered, double voxel_size, int increasing_direction, int step_direction )
{
	// Determine if on left/top or right/bottom edge, since entering a voxel can happen from either side depending on path direction, then calculate the x/y/z coordinate corresponding to the x/y/z edge, respectively
	// If stepping in direction of increasing x/y/z voxel #, entered on left/top edge, otherwise it entered on right/bottom edge.  Left/bottom edge is voxel_entered * voxel_size from 0 coordinate of first x/y/z voxel
	int on_edge = ( step_direction == increasing_direction ) ? voxel_entered : voxel_entered + 1;
	return zero_coordinate + ( increasing_direction * on_edge * voxel_size );
}
double path_projection( double m, double current_coordinate, double zero_coordinate, int current_voxel, double voxel_size, int increasing_direction, int step_direction )
{
	// Based on the dimensions of a voxel and the current (x,y,z) position, we can determine how far it is to the next edge in the x, y, and z directions.  Since the points where a path crosses 
	// one of these edges each have a corresponding (x,y,z) coordinate, we can determine which edge will be crossed next by comparing the coordinates of the next x/y/z edge in one of the three 
	// directions and determining which is closest to the current position.  For example, the x/y/z edge whose x coordinate is closest to the current x coordinate is the next edge 
	int next_voxel = current_voxel + increasing_direction * step_direction;
	double next_edge = edge_coordinate( zero_coordinate, next_voxel, voxel_size, increasing_direction, step_direction );
	// y = m(x-x0) + y0 => distance = m * (x - x0)
	return m * ( next_edge - current_coordinate );
}
double corresponding_coordinate( double m, double x, double x0, double y0 )
{
	// Using the coordinate returned by edge_coordinate, call this function to determine one of the other coordinates using 
	// y = m(x-x0)+y0 equation determine coordinates in other directions by subsequent calls to this function
	return m * ( x - x0 ) + y0;
}
/***********************************************************************************************************************************************************************************************************************/
/********************************************************************************** Proton CT Preprocessing and Image Reconstruction Code ******************************************************************************/
/***********************************************************************************************************************************************************************************************************************/
void get_dir_filenames(std::vector<std::string> &out, const std::string &directory)
{
	#if defined(_WIN32) || defined(_WIN64)
		HANDLE dir;
		WIN32_FIND_DATA file_data;
		if ((dir = FindFirstFile((directory + "/*").c_str(), &file_data)) == INVALID_HANDLE_VALUE)
			return; /* No files found */

		do {
			const std::string file_name = file_data.cFileName;
			const std::string full_file_name = directory + "/" + file_name;
			const bool is_directory = (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;

			if (file_name[0] == '.')
				continue;

			if (is_directory)
				continue;

			out.push_back(full_file_name);
		} while (FindNextFile(dir, &file_data));

		FindClose(dir);
	#else
		DIR *dir;
		class dirent *ent;
		class stat st;

		dir = opendir(directory);
		while ((ent = readdir(dir)) != NULL) {
			const string file_name = ent->d_name;
			const string full_file_name = directory + "/" + file_name;

			if (file_name[0] == '.')
				continue;

			if (stat(full_file_name.c_str(), &st) == -1)
				continue;


			const bool is_directory = (st.st_mode & S_IFDIR) != 0;

			if (is_directory)
				continue;

			out.push_back(full_file_name);
		}
		closedir(dir);
	#endif
} // get_dir_filenames
void get_dir_filenames_matching(const std::string &directory, const std::string &beginning_with, std::vector<std::string> &path, std::vector<std::string> &file)
{
	#if defined(_WIN32) || defined(_WIN64)
		HANDLE dir;
		WIN32_FIND_DATA file_data;
		if ((dir = FindFirstFile((directory + "/" + beginning_with + "*").c_str(), &file_data)) == INVALID_HANDLE_VALUE)
			return; /* No files found */

		do {
			const std::string file_name = file_data.cFileName;
			const std::string full_file_name = directory + "/" + file_name;
			const bool is_directory = (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;

			if (file_name[0] == '.')
				continue;

			if (is_directory)
				continue;

			path.push_back(full_file_name);
			file.push_back(file_name);
		} while (FindNextFile(dir, &file_data));

		FindClose(dir);
	#else
		DIR *dir;
		class dirent *ent;
		class stat st;

		dir = opendir(directory);
		while ((ent = readdir(dir)) != NULL) {
			const string file_name = ent->d_name;
			const string full_file_name = directory + "/" + file_name;

			if (file_name[0] == '.')
				continue;

			if (stat(full_file_name.c_str(), &st) == -1)
				continue;

			const bool is_directory = (st.st_mode & S_IFDIR) != 0;

			if (is_directory)
				continue;

			path.push_back(full_file_name);
			file.push_back(file_name);
		}
		closedir(dir);
	#endif
} // get_dir_filenames_matching
std::string terminal_response(char* system_command) 
{
	#if defined(_WIN32) || defined(_WIN64)
		FILE* pipe = _popen(system_command, "r");
    #else
		FILE* pipe = popen(system_command, "r");
    #endif
    
    if (!pipe) return "ERROR";
    char buffer[256];
    std::string result;
    while(!feof(pipe)) {
    	if(fgets(buffer, 256, pipe) != NULL)
    		result += buffer;
    }
	#if defined(_WIN32) || defined(_WIN64)
		 _pclose(pipe);
    #else
		 pclose(pipe);
    #endif
   
    return result;
}
void timer( bool start, clock_t start_time, clock_t end_time)
{
	if( start )
		start_time = clock();
	else
	{
		end_time = clock();
		clock_t execution_clock_cycles = (end_time - start_time) - pause_cycles;
		double execution_time = double( execution_clock_cycles) / CLOCKS_PER_SEC;
		puts("-------------------------------------------------------------------------------");
		printf( "----------------- Total execution time : %3f [seconds] -------------------\n", execution_time );	
		puts("-------------------------------------------------------------------------------");
	}
}
void exit_program_if( bool early_exit)
{
	if( early_exit )
	{
		char user_response[20];		
		timer( STOP, program_start, program_end );
		//puts("*******************************************************************************");
		//puts("************** Press enter to exit and close the console window ***************");
		//puts("*******************************************************************************");
		print_section_header( "Press 'ENTER' to exit program and close console window", '*' );
		fgets(user_response, sizeof(user_response), stdin);
		exit(1);
	}
}
void print_section_separator(char separation_char )
{
	char section_separator[CONSOLE_WINDOW_WIDTH];
	for( int n = 0; n < CONSOLE_WINDOW_WIDTH; n++ )
		section_separator[n] = separation_char;
	section_separator[CONSOLE_WINDOW_WIDTH - 1] = '\0';
	puts(section_separator);
}
void construct_header_line( char* statement, char separation_char, char* header )
{
	int length = strlen(statement), max_line_length = 70;
	int num_dashes = CONSOLE_WINDOW_WIDTH - length - 2;
	int leading_dashes = num_dashes / 2;
	int trailing_dashes = num_dashes - leading_dashes;
	int i = 0, line_length, index = 0;

	line_length = min(length - index, max_line_length);
	if( line_length < length - index )
		while( statement[index + line_length] != ' ' )
			line_length--;
	num_dashes = CONSOLE_WINDOW_WIDTH - line_length - 2;
	leading_dashes = num_dashes / 2;		
	trailing_dashes = num_dashes - leading_dashes;

	for( ; i < leading_dashes; i++ )
		header[i] = separation_char;
	if( statement[index] != ' ' )
		header[i++] = ' ';
	else
		header[i++] = separation_char;
	for( int j = 0; j < line_length; j++)
		header[i++] = statement[index++];
	header[i++] = ' ';
	for( int j = 0; j < trailing_dashes; j++)
		header[i++] = separation_char;
	header[CONSOLE_WINDOW_WIDTH - 1] = '\0';
}
void print_section_header( char* statement, char separation_char )
{
	print_section_separator(separation_char);
	char header[CONSOLE_WINDOW_WIDTH];
	int length = strlen(statement), max_line_length = 70;
	int num_dashes = CONSOLE_WINDOW_WIDTH - length - 2;
	int leading_dashes = num_dashes / 2;
	int trailing_dashes = num_dashes - leading_dashes;
	int i, line_length, index = 0;

	while( index < length )
	{
		i = 0;
		line_length = min(length - index, max_line_length);
		if( line_length < length - index )
			while( statement[index + line_length] != ' ' )
				line_length--;
		num_dashes = CONSOLE_WINDOW_WIDTH - line_length - 2;
		leading_dashes = num_dashes / 2;		
		trailing_dashes = num_dashes - leading_dashes;

		for( ; i < leading_dashes; i++ )
			header[i] = separation_char;
		if( statement[index] != ' ' )
			header[i++] = ' ';
		else
			header[i++] = separation_char;
		for( int j = 0; j < line_length; j++)
			header[i++] = statement[index++];
		header[i++] = ' ';
		for( int j = 0; j < trailing_dashes; j++)
			header[i++] = separation_char;
		header[CONSOLE_WINDOW_WIDTH - 1] = '\0';
		puts(header);		
	}
	print_section_separator(separation_char);
	puts("");
}
void print_section_exit( char* statement, char* leading_statement_chars )
{
	puts("");
	//print_section_separator('_');
	char header[CONSOLE_WINDOW_WIDTH];
	int length = strlen(statement), max_line_length = 70;
	int additional_chars = CONSOLE_WINDOW_WIDTH - length - 2;
	int leading_chars = strlen(leading_statement_chars);
	int i = 0, index = 0, line_length;

	while( index < length )
	{
		//i = 0;
		line_length = min(length - index, max_line_length);
		if( line_length < length - index )
			while( statement[index + line_length] != ' ' )
				line_length--;
		additional_chars = CONSOLE_WINDOW_WIDTH - line_length - 2;
		leading_chars = strlen(leading_statement_chars);

		if( index == 0 )
		{
			for( ; i < leading_chars; i++ )
				header[i] = leading_statement_chars[i];
		}
		else
		{
			for( ; i < leading_chars; i++ )
				header[i] = ' '; 
		}
		
		header[i++] = ' ';
		if( statement[index] == ' ' )
			index++;
		for( int j = 0; j < line_length; j++ )
			header[i++] = statement[index++];
		for( ; i < CONSOLE_WINDOW_WIDTH; i++ )
			header[i] = ' ';
		header[CONSOLE_WINDOW_WIDTH - 1] = '\0';
		puts(header);		
		i = 0;
	}
	//print_section_separator(separation_char);
	puts("");
}
char((&current_MMDD( char(&date_MMDD)[5]))[5])
{
	time_t rawtime;
	time (&rawtime);
	struct tm * timeinfo = gmtime (&rawtime);
	strftime (date_MMDD,11,"%m%d", timeinfo);
	return date_MMDD;
}
char((&current_MMDDYYYY( char(&date_MMDDYYYY)[9]))[9])
{
	time_t rawtime;
	time (&rawtime);
	struct tm * timeinfo = gmtime (&rawtime);
	strftime (date_MMDDYYYY,11,"%m%d%Y", timeinfo);
	return date_MMDDYYYY;
}
void set_IO_paths()
{
	if( true )
	{
		std::string str =  terminal_response("chdir");
		const char* cstr = str.c_str();
		current_directory = (char*) calloc( strlen(cstr), sizeof(char));
		std::copy( cstr, &cstr[strlen(cstr) - 1], current_directory );
		print_section_header( "Config file location set to current execution directory :", '*' );	
		print_section_separator('-');
		printf("%s\n", current_directory );
		print_section_separator('-');
	}
	//char* current_directory, targets_directory, WED_results_directory, target_volume_directory;
	current_MMDDYYYY( EXECUTION_DATE);
	//current_directory = (char*) calloc( strlen(cstr), sizeof(char));
	targets_directory = (char*) calloc( strlen(current_directory) + strlen(targets_folder) + 1, sizeof(char));
	WED_results_directory = (char*) calloc( strlen(current_directory) + strlen(WED_results_folder) + strlen(EXECUTION_DATE) + 1, sizeof(char));
	target_volume_directory = (char*) calloc( strlen(current_directory) + strlen(target_volume_folder) + 1, sizeof(char));

	//sprintf(current_directory, "%s\\%s", current_directory, CONFIG_FILENAME );
	sprintf(targets_directory, "%s\\%s", current_directory, targets_folder );
	sprintf(WED_results_directory, "%s\\%s_%s", current_directory, WED_results_folder, EXECUTION_DATE );
	sprintf(target_volume_directory, "%s\\%s", current_directory, target_volume_folder );
	
	char mkdir_command[256];
	sprintf(mkdir_command, "mkdir \"%s\"", WED_results_directory );
	system( mkdir_command );
	//char targets_folder[] = "bap_coordinates";
	//char target_volume_folder[] = "RStP_DICOM_PHANTOM";
	//char WED_results_folder[] = "WED_Results";
}
/***********************************************************************************************************************************************************************************************************************/
/******************************************************************************************* WED calculation and procedure functions ***********************************************************************************/
/***********************************************************************************************************************************************************************************************************************/
void WED_analysis()
{
	//combine_image_slices(target_object_dir, "_RStP_dicom_phantom_txt125mm.txt", "RStP_Phantom", WED_TARGET_ROWS, WED_TARGET_COLUMNS, WED_TARGET_SLICES, 1, BINARY);
	//combine_image_slices(target_object_dir, "_RStP_dicom_phantom_txt125mm.txt", "RStP_Phantom", WED_TARGET_ROWS, WED_TARGET_COLUMNS, WED_TARGET_SLICES, 1, TEXT);
	//read_txt_Phantom();
	read_bin_Phantom();
	find_target_beam_angles();
	for( int i = 0; i < beam_angles.size(); i++ )
	{
		read_target_coordinates( beam_angles[i] );
		WED_results = return_WED_results(beam_angles[i]);
		write_WED_results( WED_results, beam_angles[i]);
		//write_phantom_entries( BEAM_ANGLE_3 );
	}
}
void read_Phantom()
{
	float line;
	char text[20];
	char data_filename[256];
	RSP_Phantom_image_h = (float*) calloc( WED_TARGET_VOXELS, sizeof(float) );
	int total_voxels = 0;
	for( int slice = 1; slice <= WED_TARGET_SLICES; slice++ )
	{
		cout << "Slice number " << slice << endl;	
		sprintf( data_filename, "%s%03d_RStP_dicom_phantom_txt125mm.txt", target_volume_directory, slice );
		std::ifstream myfile (data_filename);
		if (myfile.is_open())
		{
			while ( myfile >> line )
				RSP_Phantom_image_h[total_voxels++] = line;	
			myfile.close();
		}			
	}
	cout << "There should be " << WED_TARGET_VOXELS/128 << " voxels in each slice" << endl;
	cout << "There should be " << WED_TARGET_VOXELS << " voxels in the phantom" << endl;
	cout << "There are " << total_voxels << " voxels in the phantom" << endl;
}
void read_txt_Phantom()
{
	char data_filename[256];
	RSP_Phantom_image_h = (float*) calloc( WED_TARGET_VOXELS, sizeof(float) );
	sprintf( data_filename, "%s\\%s.txt", target_volume_directory, target_volume_filename );
	FILE* input = fopen(data_filename, "r");
	fread(RSP_Phantom_image_h, sizeof(float), WED_TARGET_VOXELS, input );
	fclose(input);
}
void read_bin_Phantom()
{
	char data_filename[256];
	RSP_Phantom_image_h = (float*) calloc( WED_TARGET_VOXELS, sizeof(float) );
	sprintf( data_filename, "%s\\%s.bin", target_volume_directory, target_volume_filename );
	FILE* input = fopen(data_filename, "rb");
	fread(RSP_Phantom_image_h, sizeof(float), WED_TARGET_VOXELS, input );
	fclose(input);
}
//void read_slice( int slice)
//{
//	float line;
//	char text[20];
//	char data_filename[256];
//	RSP_Phantom_slice_h = (float*) calloc( WED_TARGET_COLUMNS * WED_TARGET_ROWS, sizeof(float) );
//	int num_voxels = 0;
//	sprintf( data_filename, "%s%03d_RStP_dicom_phantom_txt125mm.txt", target_object_dir, slice + 1 );
//	//cout << data_filename << endl;
//	std::ifstream myfile (data_filename);
//	if (myfile.is_open())
//	{
//		while ( myfile >> line )
//		{
//			RSP_Phantom_slice_h[num_voxels] = line;
//			num_voxels++;		
//		}
//		myfile.close();
//	}
//	cout << "Slice number " << slice << endl;
//	cout << "Elements read = " << num_voxels << endl;
//	cout << "There should be " << WED_TARGET_VOXELS/128 << " voxels in this slice" << endl;
//	//fgets(text, sizeof text, stdin);
//	//write_slice( "RSP_Phantom2", target_object_dir, RSP_Phantom_slice_h, WED_TARGET_COLUMNS, WED_TARGET_ROWS, slice);
//}
void write_Phantom()
{
	char output_filename_txt[256];
	sprintf(output_filename_txt, "%s\\%s.txt", target_volume_directory, target_volume_filename);
	FILE* combined_image_file_txt = fopen(output_filename_txt, "w");
	fwrite(RSP_Phantom_image_h, sizeof(float), WED_TARGET_VOXELS, combined_image_file_txt);
	fclose(combined_image_file_txt);

	char output_filename_bin[256];
	sprintf(output_filename_bin, "%s\\%s.bin", target_volume_directory, target_volume_filename);
	FILE* combined_image_file_bin = fopen(output_filename_bin, "wb");
	fwrite(RSP_Phantom_image_h, sizeof(float), WED_TARGET_VOXELS, combined_image_file_bin);
	fclose(combined_image_file_bin);
}
void write_slice( char* output_filename_base, const char* output_directory, float* &float_array, int x_max, int y_max, int slice )
{
	char output_filename[256];
	// Write each slice of the array/image to a separate file
	std::ofstream output_file;
	sprintf( output_filename, "%s%s_%d.txt", output_directory, output_filename_base, slice );
	output_file.open(output_filename);
	for(int y = 0; y < y_max; y++)
	{
		for(int x = 0; x < x_max; x++)
			output_file << float_array[ ( y * x_max ) + x ] << " ";
		output_file << endl;
	}
	output_file.close();	
}
void combine_image_slices(const char* directory, char* file_basename, char* combined_filename, uint rows, uint columns, uint slices, uint start_slice, DISK_WRITE_MODE format )
{
	float line;
	char text[20];
	char data_filename[256];
	uint num_voxels = rows * columns * slices;
	float* image = (float*) calloc( num_voxels, sizeof(float) );
	int elements = 0;
	int total_voxels = 0;
	for( int slice = start_slice; slice < start_slice + slices; slice++ )
	{
		sprintf( data_filename, "%s%03d%s", directory, slice, file_basename );
		std::ifstream slice_file (data_filename);
		if (slice_file.is_open())
		{
			while ( slice_file >> line )
				image[total_voxels++] = line;
			slice_file.close();
		}
		cout << "Finished reading slice " << slice << endl;
	}
	char output_filename_txt[256];
	//if( format == TEXT )
	//{
		sprintf(output_filename_txt, "%s\\%s.txt", directory, combined_filename);
		FILE* combined_image_file_txt = fopen(output_filename_txt, "w");
		fwrite(image, sizeof(float), num_voxels, combined_image_file_txt);
		fclose(combined_image_file_txt);
	//}
		char output_filename_bin[256];
	//else if( format == BINARY )
	//{
		sprintf(output_filename_bin, "%s\\%s.bin", directory, combined_filename);
		FILE* combined_image_file_bin = fopen(output_filename_bin, "wb");
		fwrite(image, sizeof(float), num_voxels, combined_image_file_bin);
		fclose(combined_image_file_bin);
	//}
} 
void find_target_beam_angles()
{
	std::vector<std::string> path;
	std::vector<std::string> file;
	get_dir_filenames_matching(std::string(targets_directory), std::string(targets_base_name), path, file);
	int beam_angle_int;
	char * cstr = (char*)calloc(sizeof(char), file[0].length()+1);
	char * beam_angle = (char*)calloc(sizeof(char), 4);
	const char* file_name;
	char* h1, *h2;
	for( int i = 0; i < file.size(); i++ )
		cout << file[i] << endl;
	for( int i = 0; i < file.size(); i++ )
	{		
		std::strcpy (cstr, file[i].c_str());
		printf("%s\n", cstr );
		file_name = file[i].c_str();
		h1 = strstr( cstr, "_" )   + 1;
		h2 = strstr( h1, "_" )     + 1;	
		strncpy ( beam_angle, h2, 3 );
		beam_angle_int = stoi(std::string(beam_angle));
		cout << "beam angle = " << beam_angle_int  << endl;
		beam_angles.push_back(beam_angle_int);
	}
}
void read_target_coordinates( int angle )
{
	double coordinate;
	char text[20];
	char data_filename[256];
	std::vector<double> coordinates;
	sprintf( data_filename, "%s\\%s_%03d_deg.txt", targets_directory, targets_base_name, angle );
	std::ifstream myfile (data_filename);
	if (myfile.is_open())
	{
		while ( myfile >> coordinate )
			coordinates.push_back(coordinate);
		myfile.close();
	}
	num_targets = coordinates.size()/3;
	target_x_h = (double*) calloc( num_targets, sizeof(double) );
	target_y_h = (double*) calloc( num_targets, sizeof(double) );
	target_z_h = (double*) calloc( num_targets, sizeof(double) );
	for( int i = 0, j = 0; j < num_targets; i += 3, j++ )
	{
		target_x_h[j] = coordinates[i];
		target_y_h[j] = coordinates[i+1];
		target_z_h[j] = coordinates[i+2];

	}
	//print_targets();
}
void print_targets()
{
	for( int i = 0; i < num_targets; i++ )
		cout << "Target " << i << ":" << target_x_h[i] << " " << target_y_h[i] << " " << target_z_h[i] << " " << endl;
}
//double*& return_WED_results_old( int beam_angle)
//{
//	double* WED_values = (double*) calloc( num_targets, sizeof(double) );
//	int slice;
//	for( int target = 0; target < num_targets; target++ )
//	{
//		cout << "processing target # " << target + 1 << " of " << num_targets << " for beam angle " << beam_angle << endl;;
//		if( !FULL_PHANTOM )
//		{
//			slice = calculate_slice( target_z_h[target] );
//			cout << "Slice read = " << slice << endl;
//			read_slice( slice );
//		}
//		WED_values[target] = calculate_target_WED_old( target_x_h[target], target_y_h[target], target_z_h[target], beam_angle );
//	}
//	return WED_values;
//}
double*& return_WED_results( int beam_angle)
{
	double* WED_values = (double*) calloc( num_targets, sizeof(double) );
	for( int target = 0; target < num_targets; target++ )
	{
		cout << "processing target # " << target + 1 << " of " << num_targets << " for beam angle " << beam_angle << endl;;
		WED_values[target] =  calculate_target_WED( target_x_h[target], target_y_h[target], target_z_h[target], beam_angle, 0.0 );
	}
	return WED_values;
}
void WED_take_2D_step
( 
	const int x_move_direction, const int y_move_direction, const int z_move_direction, 
	const double dy_dx, const double dz_dx, const double dz_dy, const double dx_dy, const double dx_dz, const double dy_dz, 
	const double x_start, const double y_start, const double z_start, double& x, double& y, double& z, 
	int& voxel_x, int& voxel_y, int& voxel_z, int& voxel, double& x_to_go, double& y_to_go, double& z_to_go, double& chord_length	
)
{
	double previous_x = x, previous_y = y, previous_z = z;
	// Change in x for Move to Voxel Edge in y
	double y_extension = fabs( dx_dy ) * y_to_go;
	//If Next Voxel Edge is in x or xy Diagonal
	if( x_to_go <= y_extension )
	{
		//printf(" x_to_go <= y_extension \n");
		voxel_x += x_move_direction;					
		x = edge_coordinate( WED_TARGET_X_ZERO_COORDINATE, voxel_x, WED_TARGET_VOXEL_WIDTH, X_INCREASING_DIRECTION, x_move_direction );
		y = corresponding_coordinate( dy_dx, x, x_start, y_start );
		x_to_go = WED_TARGET_VOXEL_WIDTH;
		y_to_go = distance_remaining( WED_TARGET_Y_ZERO_COORDINATE, y, Z_INCREASING_DIRECTION, y_move_direction, WED_TARGET_VOXEL_HEIGHT, voxel_y );
	}
	// Else Next Voxel Edge is in y
	else
	{
		//printf(" y_extension < x_extension \n");				
		voxel_y -= y_move_direction;
		y = edge_coordinate( WED_TARGET_Y_ZERO_COORDINATE, voxel_y, WED_TARGET_VOXEL_HEIGHT, Y_INCREASING_DIRECTION, y_move_direction );
		x = corresponding_coordinate( dx_dy, y, y_start, x_start );
		x_to_go = distance_remaining( WED_TARGET_X_ZERO_COORDINATE, x, X_INCREASING_DIRECTION, x_move_direction, WED_TARGET_VOXEL_WIDTH, voxel_x );
		y_to_go = WED_TARGET_VOXEL_HEIGHT;
	}
	if( x_to_go == 0 )
	{
		x_to_go = WED_TARGET_VOXEL_WIDTH;
		voxel_x += x_move_direction;
	}
	if( y_to_go == 0 )
	{
		y_to_go = WED_TARGET_VOXEL_HEIGHT;
		voxel_y -= y_move_direction;
	}
	voxel_z = max(voxel_z, 0 );
	voxel = voxel_x + voxel_y * WED_TARGET_COLUMNS + voxel_z * WED_TARGET_COLUMNS * WED_TARGET_ROWS;
	double delta_x_sqd = pow(x - previous_x, 2.0);
	double delta_y_sqd = pow(y - previous_y, 2.0);
	double delta_z_sqd = pow(z - previous_z, 2.0);
	chord_length = sqrt( delta_x_sqd + delta_y_sqd + delta_z_sqd );
}
void WED_take_3D_step
( 
	const int x_move_direction, const int y_move_direction, const int z_move_direction,
	const double dy_dx, const double dz_dx, const double dz_dy, const double dx_dy, const double dx_dz, const double dy_dz, 
	const double x_start, const double y_start, const double z_start, double& x, double& y, double& z, 
	int& voxel_x, int& voxel_y, int& voxel_z, int& voxel, double& x_to_go, double& y_to_go, double& z_to_go, double& chord_length		
)
{
	double previous_x = x, previous_y = y, previous_z = z;
	// Change in z for Move to Voxel Edge in x and y
	double x_extension = fabs( dz_dx ) * x_to_go;
	double y_extension = fabs( dz_dy ) * y_to_go;
	if( (z_to_go <= x_extension  ) && (z_to_go <= y_extension) )
	{
		//printf("z_to_go <= x_extension && z_to_go <= y_extension\n");				
		voxel_z -= z_move_direction;					
		z = edge_coordinate( WED_TARGET_Z_ZERO_COORDINATE, voxel_z, WED_TARGET_VOXEL_THICKNESS, Z_INCREASING_DIRECTION, z_move_direction );					
		x = corresponding_coordinate( dx_dz, z, z_start, x_start );
		y = corresponding_coordinate( dy_dz, z, z_start, y_start );
		x_to_go = distance_remaining( WED_TARGET_X_ZERO_COORDINATE, x, X_INCREASING_DIRECTION, x_move_direction, WED_TARGET_VOXEL_WIDTH, voxel_x );
		y_to_go = distance_remaining( WED_TARGET_Y_ZERO_COORDINATE, y, Y_INCREASING_DIRECTION, y_move_direction, WED_TARGET_VOXEL_HEIGHT, voxel_y );	
		z_to_go = WED_TARGET_VOXEL_THICKNESS;
	}
	//If Next Voxel Edge is in x or xy Diagonal
	else if( x_extension <= y_extension )
	{
		//printf(" x_extension <= y_extension \n");					
		voxel_x += x_move_direction;
		x = edge_coordinate( WED_TARGET_X_ZERO_COORDINATE, voxel_x, WED_TARGET_VOXEL_WIDTH, X_INCREASING_DIRECTION, x_move_direction );
		y = corresponding_coordinate( dy_dx, x, x_start, y_start );
		z = corresponding_coordinate( dz_dx, x, x_start, z_start );
		x_to_go = WED_TARGET_VOXEL_WIDTH;
		y_to_go = distance_remaining( WED_TARGET_Y_ZERO_COORDINATE, y, Y_INCREASING_DIRECTION, y_move_direction, WED_TARGET_VOXEL_HEIGHT, voxel_y );
		z_to_go = distance_remaining( WED_TARGET_Z_ZERO_COORDINATE, z, Z_INCREASING_DIRECTION, z_move_direction, WED_TARGET_VOXEL_THICKNESS, voxel_z );
	}
	// Else Next Voxel Edge is in y
	else
	{
		//printf(" y_extension < x_extension \n");
		voxel_y -= y_move_direction;					
		y = edge_coordinate( WED_TARGET_Y_ZERO_COORDINATE, voxel_y, WED_TARGET_VOXEL_HEIGHT, Y_INCREASING_DIRECTION, y_move_direction );
		x = corresponding_coordinate( dx_dy, y, y_start, x_start );
		z = corresponding_coordinate( dz_dy, y, y_start, z_start );
		x_to_go = distance_remaining( WED_TARGET_X_ZERO_COORDINATE, x, X_INCREASING_DIRECTION, x_move_direction, WED_TARGET_VOXEL_WIDTH, voxel_x );
		y_to_go = WED_TARGET_VOXEL_HEIGHT;					
		z_to_go = distance_remaining( WED_TARGET_Z_ZERO_COORDINATE, z, Z_INCREASING_DIRECTION, z_move_direction, WED_TARGET_VOXEL_THICKNESS, voxel_z );
	}
	if( x_to_go == 0 )
	{
		x_to_go = WED_TARGET_VOXEL_WIDTH;
		voxel_x += x_move_direction;
	}
	if( y_to_go == 0 )
	{
		y_to_go = WED_TARGET_VOXEL_HEIGHT;
		voxel_y -= y_move_direction;
	}
	if( z_to_go == 0 )
	{
		z_to_go = WED_TARGET_VOXEL_THICKNESS;
		voxel_z -= z_move_direction;
	}
	voxel_z = max(voxel_z, 0 );
	voxel = voxel_x + voxel_y * WED_TARGET_COLUMNS + voxel_z * WED_TARGET_COLUMNS * WED_TARGET_ROWS;
	double delta_x_sqd = pow(x - previous_x, 2.0);
	double delta_y_sqd = pow(y - previous_y, 2.0);
	double delta_z_sqd = pow(z - previous_z, 2.0);
	chord_length = sqrt( delta_x_sqd + delta_y_sqd + delta_z_sqd );
}
double calculate_target_WED( double x_target, double y_target, double z_target, double beam_angle_XY, double beam_angle_XZ )
{
		char text[20];
		/****************************************************************************************************************************************/
		/****************************************************** Path Characteristic Parameters **************************************************/
		/****************************************************************************************************************************************/
		double beam_angle_XY_radians = beam_angle_XY * ANGLE_TO_RADIANS;
		double beam_angle_XZ_radians = beam_angle_XZ * ANGLE_TO_RADIANS;
		double dy_dx, dz_dx, dz_dy, dx_dy, dx_dz, dy_dz;
		double x, y, z;	
		double x_to_go, y_to_go, z_to_go;	
		int x_move_direction, y_move_direction, z_move_direction;
		int voxel_x, voxel_y, voxel_z, voxel;
		int voxel_x_out, voxel_y_out, voxel_z_out, voxel_out; 
		/****************************************************************************************************************************************/
		/******************************************************* Status Tracking Information ****************************************************/
		/****************************************************************************************************************************************/
		double RSP;
		double chord_length;
		double WED = 0;						
		bool end_walk;	
		/****************************************************************************************************************************************/
		/********************************************** Initial Conditions and Movement Characteristics *****************************************/
		/****************************************************************************************************************************************/
		double x_start, y_start, z_start;
		bool x_start_direction = cos( beam_angle_XY_radians) >= 0;
		double x_intercept_1 , x_intercept_2, y_intercept_1, y_intercept_2;
		bool x_diff;

		double m_in = tan( beam_angle_XY_radians );	// proton entry path slope
		double b_in = y_target - m_in * x_target;				// proton entry path y-intercept
		
		// Quadratic formula coefficients
		double a = 1 + pow(m_in, 2);							// x^2 coefficient 
		double b = 2 * m_in * b_in;							// x coefficient
		double c = pow(b_in, 2) - pow(WED_TARGET_IMAGE_WIDTH/2, 2.0 );		// 1 coefficient
		double entry_discriminant = pow(b, 2) - (4 * a * c);	// Quadratic formula discriminant		
		bool entered = ( entry_discriminant > 0 );			// Proton path intersected twice
		
		// Find both intersection points of the circle; closest one to the entry SSDs is the entry position
		if( entered )
		{
			// -b +- sqrt(b^2 - 4ac) => x1 = -(sqrt(b^2 - 4ac) + b)
			//							x2 = sqrt(b^2 - 4ac) -b 
			//y1 = m*x1+bin =  -m*(sqrt(b^2 - 4ac) + b) + bin
			//y2 = m*x2 + bin = 
			x_intercept_1 =   ( sqrtf(entry_discriminant) - b ) / ( 2 * a );
			x_intercept_2 = - ( sqrtf(entry_discriminant) + b ) / ( 2 * a );
			y_intercept_1 = m_in * x_intercept_1 + b_in;
			y_intercept_2 = m_in * x_intercept_2 + b_in;

			x_diff = x_intercept_1 > x_intercept_2;
			x_start = ( (cos( beam_angle_XY_radians) >= 0) &&  (x_intercept_1 > x_intercept_2)) ? x_intercept_1 : x_intercept_2;
			y_start = (x_start_direction == x_diff) ? y_intercept_1 : y_intercept_2;

			z_start = z_target;
		}	
		//cout << "x_start = " << x_start << " y_start = " << y_start << " z_start = " << z_start << endl;
		/****************************************************************************************************************************************/
		/******************************************* Initial Conditions and Movement Characteristics ********************************************/
		/****************************************************************************************************************************************/
		x = x_start;
		y = y_start;
		z = z_start;

		x_move_direction = ( x_target >= x_start ) - ( x_target <= x_start );
		y_move_direction = ( y_target >= y_start ) - ( y_target <= y_start );
		z_move_direction = ( z_target >= z_start ) - ( z_target <= z_start );

		voxel_x = calculate_voxel( WED_TARGET_X_ZERO_COORDINATE, x, WED_TARGET_VOXEL_WIDTH );
		voxel_y = calculate_voxel( WED_TARGET_Y_ZERO_COORDINATE, y, WED_TARGET_VOXEL_HEIGHT );
		voxel_z = calculate_voxel( WED_TARGET_Z_ZERO_COORDINATE, z, WED_TARGET_VOXEL_THICKNESS );
		voxel = int(voxel_x + voxel_y * WED_TARGET_COLUMNS + voxel_z * WED_TARGET_COLUMNS * WED_TARGET_ROWS);	
		cout << "voxel_x = " << voxel_x << " voxel_y = " << voxel_y << " voxel_z = " << voxel_z << endl;
		
		x_to_go = distance_remaining( WED_TARGET_X_ZERO_COORDINATE, x, X_INCREASING_DIRECTION, x_move_direction, WED_TARGET_VOXEL_WIDTH, voxel_x );
		y_to_go = distance_remaining( WED_TARGET_Y_ZERO_COORDINATE, y, Y_INCREASING_DIRECTION, y_move_direction, WED_TARGET_VOXEL_HEIGHT, voxel_y );	
		z_to_go = distance_remaining( WED_TARGET_Z_ZERO_COORDINATE, z, Z_INCREASING_DIRECTION, z_move_direction, WED_TARGET_VOXEL_THICKNESS, voxel_z );
		
		voxel_x_out = calculate_voxel( WED_TARGET_X_ZERO_COORDINATE, x_target, WED_TARGET_VOXEL_WIDTH );
		voxel_y_out = calculate_voxel( WED_TARGET_Y_ZERO_COORDINATE, y_target, WED_TARGET_VOXEL_HEIGHT );
		voxel_z_out = calculate_voxel( WED_TARGET_Z_ZERO_COORDINATE, z_target, WED_TARGET_VOXEL_THICKNESS );			
		voxel_out = int(voxel_x_out + voxel_y_out * WED_TARGET_COLUMNS + voxel_z_out * WED_TARGET_COLUMNS * WED_TARGET_ROWS);
		cout << "voxel_x_out = " << voxel_x_out << " voxel_y_out = " << voxel_y_out << " voxel_z_out = " << voxel_z_out << endl;

		end_walk = ( voxel == voxel_out ) || ( voxel_x >= WED_TARGET_COLUMNS ) || ( voxel_y >= WED_TARGET_ROWS ) || ( voxel_z >= WED_TARGET_SLICES );
		/****************************************************************************************************************************************/
		/****************************************************** Path and Walk Information *******************************************************/
		/****************************************************************************************************************************************/
		// Slopes corresponging to each possible direction/reference.  Explicitly calculated inverses to avoid 1/# calculations later
		dy_dx = tan(beam_angle_XY_radians);
		dz_dx = tan(beam_angle_XZ_radians);
		dz_dy = tan(beam_angle_XZ_radians)/tan(beam_angle_XY_radians);
		dx_dy = pow(tan(beam_angle_XY_radians), -1.0);
		dx_dz = pow(tan(beam_angle_XZ_radians), -1.0);
		dy_dz = tan(beam_angle_XY_radians)/tan(beam_angle_XZ_radians);	
		/****************************************************************************************************************************************/
		/************************************************************ Voxel Walk Routine ********************************************************/
		/****************************************************************************************************************************************/
		int i = 0;
		if( z_move_direction != 0 )
		{
			//printf("z_exit[i] != z_entry[i]\n");
			while( !end_walk )
			{
				RSP = RSP_Phantom_image_h[voxel];
				WED_take_3D_step
				( 
					x_move_direction, y_move_direction, z_move_direction, dy_dx, dz_dx, dz_dy, dx_dy, dx_dz, dy_dz, 
					x_start, y_start, z_start, x, y, z, voxel_x, voxel_y, voxel_z, voxel, x_to_go, y_to_go, z_to_go, chord_length
				);	
				if( RSP >= WED_TARGET_THRESHOLD_RSP )
					WED += chord_length * RSP;
				end_walk = ( voxel == voxel_out ) || ( voxel_x >= WED_TARGET_COLUMNS ) || ( voxel_y >= WED_TARGET_ROWS ) || ( voxel_z >= WED_TARGET_SLICES );
				i++;
			}// end !end_walk 
		}
		else
		{
			//printf("z_exit[i] == z_entry[i]\n");
			while( !end_walk )
			{
				RSP = RSP_Phantom_image_h[voxel];
				WED_take_2D_step
				( 
					x_move_direction, y_move_direction, z_move_direction, dy_dx, dz_dx, dz_dy, dx_dy, dx_dz, dy_dz, x_start, y_start, z_start, 
					x, y, z, voxel_x, voxel_y, voxel_z, voxel, x_to_go, y_to_go, z_to_go, chord_length
				);							
				if( RSP >= WED_TARGET_THRESHOLD_RSP )
					WED += chord_length * RSP;
				end_walk = ( voxel == voxel_out ) || ( voxel_x >= WED_TARGET_COLUMNS ) || ( voxel_y >= WED_TARGET_ROWS ) || ( voxel_z >= WED_TARGET_SLICES );
				i++;
			}// end: while( !end_walk )
		}//end: else: z_start != z_end => z_start == z_end
		/****************************************************************************************************************************************/
		/************************************************ Step from final voxel edge to target **************************************************/
		/****************************************************************************************************************************************/
		//cout << "voxels passed through = " << i << endl;
		//cout << "WED before last partial step = " << WED << endl;
		double delta_x_sqd = pow(x - x_target, 2.0);
		double delta_y_sqd = pow(y - y_target, 2.0);
		double delta_z_sqd = pow(z - z_target, 2.0);
		chord_length = sqrt( delta_x_sqd + delta_y_sqd + delta_z_sqd );
		WED += chord_length* RSP_Phantom_image_h[voxel];
		//cout << "voxels passed through at least partially = " << i + 1<< endl;
		cout << "WED after final partial step = " << WED << endl;	
		return WED;
}
void write_phantom_entries(int beam_angle)
{
	char output_filename[256];
	std::ofstream output_file;
	sprintf( output_filename, "%s\\Entries_%03d.txt", WED_results_directory, beam_angle );
	output_file.open(output_filename);
	for( int i = 0; i < entry_voxel_x.size(); i++ )
		output_file << entry_voxel_x[i] << " " <<  entry_voxel_y[i] << " " <<  entry_voxel_z[i] << endl;
	output_file.close();
	cout << entry_voxel_x.size() << endl;
}
void write_WED_results( double*& WED_values, int beam_angle)
{
	char output_filename[256];
	std::ofstream output_file;
	sprintf( output_filename, "%s\\%s_%03d.txt", WED_results_directory, WED_results_basename, beam_angle );
	output_file.open(output_filename);
	//cout << "num_targets = " << num_targets << endl;
	for(int target = 0; target < num_targets; target++)
	{
		output_file << WED_values[target] << endl;
	}
	output_file.close();	
}