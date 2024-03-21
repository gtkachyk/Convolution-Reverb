#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <arpa/inet.h>
#include <algorithm>
#include <climits>
#include <iostream>
#include <math.h>

// Constants
#define MASK 0xFF
#define BITS_PER_BYTE 8
#define STRING_FIELD_SIZE 4
#define SUBCHUNK_1_SIZE_NORMAL_VALUE 16
#define CHUNK_1_SIZE 36
#define PCM_AUDIO_FORMAT_IDENTIFIER 1
#define BITS_PER_SAMPLE_DEFAULT 16
#define CHUNK1ID_DEFAULT "RIFF"
#define FILE_TYPE_HEADER_DEFAULT "WAVE"
#define SUBCHUNK_1_ID_DEFAULT "fmt "
#define SUBCHUNK_2_ID_DEFAULT "data"

// FFT Data.
#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr
#define FFT_ISIGN 1
#define INVERSE_FFT_ISIGN -1

// Macros.
#define INT_TO_BIG_ENDIAN(int_to_format, current_byte) (int_to_format >> (24 - BITS_PER_BYTE * current_byte)) & MASK
#define INT_TO_LITTLE_ENDIAN(int_to_format, current_byte) (int_to_format >> (BITS_PER_BYTE * current_byte)) & MASK
#define SHORT_TO_BIG_ENDIAN(short_to_format, current_byte) (short_to_format >> (BITS_PER_BYTE - BITS_PER_BYTE * current_byte)) & MASK
#define SHORT_TO_LITTLE_ENDIAN(short_to_format, current_byte) (short_to_format >> (BITS_PER_BYTE * current_byte)) & MASK

bool is_big_endian;

struct wav_file {
    // Header:
    // RIFF chunk descriptor.
    char chunk_id[STRING_FIELD_SIZE];
    uint32_t chunk_size;
    char format[STRING_FIELD_SIZE];

    // fmt sub-chunk.
    char subchunk_1_id[STRING_FIELD_SIZE];
    uint32_t subchunk_1_size;
    uint16_t audio_format;
    uint16_t num_channels;
    uint32_t sample_rate;
    uint32_t byte_rate;
    uint16_t block_align;
    uint16_t bits_per_sample;

    // data sub-chunk.
    char subchunk_2_id[STRING_FIELD_SIZE];
    uint32_t subchunk_2_size;
    char* data;

    // Additional data.
    short* signal;
    int signal_size;
};

/**
 * Prints the fields of a wav_file struct.
 * @param wav_file The wav_file struct to print.
 * @param wav_file_name The name of 'wav_file.'
 */
void print_wav_file(wav_file wav_file, std::string wav_file_name){
    // RIFF chunk descriptor.
    printf("%s.chunk_id = %s\n", wav_file_name.c_str(), wav_file.chunk_id);
    printf("%s.chunk_size = %d\n", wav_file_name.c_str(), wav_file.chunk_size);
    printf("%s.format = %s\n", wav_file_name.c_str(), wav_file.format);

    // fmt sub-chunk.
    printf("%s.subchunk_1_id = %s\n", wav_file_name.c_str(), wav_file.subchunk_1_id);
    printf("%s.subchunk_1_size = %d\n", wav_file_name.c_str(), wav_file.subchunk_1_size);
    printf("%s.audio_format = %d\n", wav_file_name.c_str(), wav_file.audio_format);
    printf("%s.num_channels = %d\n", wav_file_name.c_str(), wav_file.num_channels);
    printf("%s.sample_rate = %d\n", wav_file_name.c_str(), wav_file.sample_rate);
    printf("%s.byte_rate = %d\n", wav_file_name.c_str(), wav_file.byte_rate);
    printf("%s.block_align = %d\n", wav_file_name.c_str(), wav_file.block_align);
    printf("%s.bits_per_sample = %d\n", wav_file_name.c_str(), wav_file.bits_per_sample);

    // data sub-chunk.
    printf("%s.subchunk_2_id = %s\n", wav_file_name.c_str(), wav_file.subchunk_2_id);
    printf("%s.subchunk_2_size = %d\n", wav_file_name.c_str(), wav_file.subchunk_2_size);
    printf("%s.signal_size = %d\n", wav_file_name.c_str(), wav_file.signal_size);

}

/**
 * Reads data of a given type from the provided file stream.
 * @tparam T The data type of the value to be read from the file.
 * @param file An input file stream representing the opened file.
 * @param value A reference to the variable where the read value will be stored.
 */
template <typename T> void read_from_stream(std::ifstream& file, T& value) {
    // Read data from the file stream and store it in the specified variable.
    file.read(reinterpret_cast<char*>(&value), sizeof(T));
}

/**
 * Reads any extra bytes if the sub-chunk size in the 'fmt' sub-chunk is equal to 18.
 * @param file An input file stream representing the opened WAV file.
 * @param subchunk_1_size The size of the 'fmt' sub-chunk in the WAV file.
 */
void read_extra_bytes(std::ifstream& file, uint32_t subchunk_1_size) {
    if (subchunk_1_size == SUBCHUNK_1_SIZE_NORMAL_VALUE + 2) {
        // If extra bytes are present, read and discard them.
        short empty_bytes;
        read_from_stream(file, empty_bytes);
    }
}

/**
 * Reads the data sub-chunk of a WAV file from the provided file stream.
 * @param file An input file stream representing the opened WAV file.
 * @param wav_file_struct The WAV file structure where the read data will be stored.
 */
void read_data_subchunk(std::ifstream& file, wav_file& wav_file_struct) {
    // Read sub-chunk ID and sub-chunk size.
    read_from_stream(file, wav_file_struct.subchunk_2_id);
    read_from_stream(file, wav_file_struct.subchunk_2_size);

    // Calculate the size of the audio data.
    int data_size = wav_file_struct.subchunk_2_size;

    // Allocate memory for the raw audio data.
    wav_file_struct.data = new char[data_size];

    // Read the raw audio data from the file.
    file.read(wav_file_struct.data, data_size);
}

/**
 * Converts the raw audio data in a WAV file structure to a short signal based on the bits per sample.
 * Modified from https://stackoverflow.com/questions/66225400/how-to-merge-two-unsigned-chars-into-a-single-unsigned-short16-bits-in-c.
 * @param wav_file_struct The WAV file structure containing the raw audio data.
 */
void convert_signal(wav_file& wav_file_struct) {
    // For 16-bit audio, pairs of consecutive unsigned chars are combined into a short.
    if(wav_file_struct.bits_per_sample == BITS_PER_SAMPLE_DEFAULT){
        int dataSize = wav_file_struct.subchunk_2_size;
        wav_file_struct.signal_size = dataSize / 2;
        wav_file_struct.signal = new short[wav_file_struct.signal_size];
        short signal_element;
        for(int i = 0; i < dataSize; i+=2){
            signal_element = static_cast<short>(static_cast<unsigned char> (wav_file_struct.data[i]));
            signal_element = static_cast<short>(static_cast<unsigned char> (wav_file_struct.data[i + 1])) << BITS_PER_BYTE;
            wav_file_struct.signal[i / 2] = signal_element;
        }
    }
    else{
        printf("Error: input wave files must have a bit depth of 16. Exiting...\n");
        exit(1);
    }
}

/**
 * Creates a new wav_file struct by reading the contents of a WAV file.
 * @param file_name The name of the WAV file to be read.
 * @return A 'wav_file' structure containing the metadata and converted audio signal of the input WAV file.
 */
wav_file read_wav_file(const std::string& file_name) {
    // Create a new wav_file struct to store the WAV file data.
    wav_file wav_file_struct;

    // Open the WAV file.
    std::ifstream file(file_name, std::ios::in | std::ios::binary);

    // Read the header fields from the file.
    read_from_stream(file, wav_file_struct.chunk_id);
    read_from_stream(file, wav_file_struct.chunk_size);
    read_from_stream(file, wav_file_struct.format);

    read_from_stream(file, wav_file_struct.subchunk_1_id);
    read_from_stream(file, wav_file_struct.subchunk_1_size);
    read_from_stream(file, wav_file_struct.audio_format);
    read_from_stream(file, wav_file_struct.num_channels);
    read_from_stream(file, wav_file_struct.sample_rate);
    read_from_stream(file, wav_file_struct.byte_rate);
    read_from_stream(file, wav_file_struct.block_align);
    read_from_stream(file, wav_file_struct.bits_per_sample);

    // Handle any extra bytes in the 'fmt' sub-chunk.
    read_extra_bytes(file, wav_file_struct.subchunk_1_size);

    // Read the data sub-chunk and store the raw audio data.
    read_data_subchunk(file, wav_file_struct);

    // Convert the raw audio data into a short signal based on the bits per sample.
    convert_signal(wav_file_struct);

    return wav_file_struct;
}

/**
 * Normalize a set of short integer values to the range [-1.0, 1.0].
 * @param original_signal The signal to prepare.
 * @param signal_size The size of the signal to prepare.
 */
void prepare_signal_for_convolving(const short* original_signal, double new_signal[], int signal_size) {
    for (int i = 0; i < signal_size; i++) {
        new_signal[i] = static_cast<double>(original_signal[i]) / SHRT_MAX;
    }
}

/**
 * Modifies the values in a convolved signal to match the magnitude of the original signal.
 * @param convolved_signal Pointer to the array containing the output signal to be scaled.
 * @param original_file Pointer to a wav_file struct representing the input WAV file, containing the original signal to determine the scaling factor.
 * @param new_signal_size The size of the output signal array.
 */
void set_convolved_signal_magnitude(wav_file* original_file, double* convolved_signal, int new_signal_size) {
    double input_max = *std::max_element((*original_file).signal, (*original_file).signal + new_signal_size);
    double output_max = *std::max_element(convolved_signal, convolved_signal + new_signal_size);

    // Avoid division by zero and handle potential precision issues.
    if (output_max != 0.0) {
        double scaler = input_max / output_max;

        for (int i = 0; i < new_signal_size; ++i) {
            convolved_signal[i] *= scaler;
        }
    } 
    else {
        // Handle the case where output_max is zero.
        std::fill(convolved_signal, convolved_signal + new_signal_size, 0.0);
    }
}

/**
 * Formats an int into big endian or little endian depending on the value of 'is_big_endian.'
 * @param int_to_format The int to format.
 * @param formatted_int An array of unsigned characters used to store the formatted int.
 */
void format_int_endianness(int int_to_format, unsigned char* formatted_int) {
    int size_of_int = (int)sizeof(int);
    if(is_big_endian){
        for (int i = 0; i < size_of_int; i++) {
            formatted_int[i] = static_cast<unsigned char>(INT_TO_BIG_ENDIAN(int_to_format, i));
        }
    }
    else{
        for (int i = 0; i < size_of_int; i++) {
            formatted_int[i] = static_cast<unsigned char>(INT_TO_LITTLE_ENDIAN(int_to_format, i));
        }
    }
}

/**
 * Formats a short into big endian or little endian depending on the value of 'is_big_endian.'
 * @param short_to_format The short to format.
 * @param formatted_short An array of unsigned characters used to store the formatted short.
 */
void format_short_endianness(short short_to_format, unsigned char* formatted_short) {
    int size_of_short = (int)sizeof(short);
    if(is_big_endian){
        for (int i = 0; i < size_of_short; i++) {
            formatted_short[i] = static_cast<unsigned char>(SHORT_TO_BIG_ENDIAN(short_to_format, i));
        }
    }
    else{
        for (int i = 0; i < size_of_short; i++) {
            formatted_short[i] = static_cast<unsigned char>(SHORT_TO_LITTLE_ENDIAN(short_to_format, i));
        }
    }
}

/**
 * Sets the value of 'is_big_endian' global variable.
 * Taken from https://stackoverflow.com/questions/1001307/detecting-endianness-programmatically-in-a-c-program.
 */
void set_endianness(){
    int test_int = 47;
    if(htonl(test_int) == test_int) {
        is_big_endian = true;
    } 
    else {
        is_big_endian = false;
    }
}

/**
 * The four1 FFT from Numerical Recipes in C, p. 507 - 508. 
 * Note: changed float data types to double.
 * Parameter nn must be a power of 2, and use +1 for isign for an FFT, and -1 for the Inverse FFT.
 * The data is complex, so the array size must be nn*2. 
 * This code assumes the array starts at index 1, not 0, so subtract 1 when calling the routine (see main() below).
 * @param data A complex array representing the input and output data. The array must have a size of nn*2, and the indexing starts at 1.
 * @param nn The size of the data array, which must be a power of 2.
 * @param isign The sign of the transform. A value of 1 is used for FFT and -1 is used for Inverse FFT.
 */
void four1(double data[], int nn, int isign){
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}

/**
 * Finds the largest power of 2 that is less than or equal to a specified value.
 * Taken from https://www.geeksforgeeks.org/highest-power-2-less-equal-given-number/.
 * @param n The limiting value.
 * @return The largest power of 2 that is less than or equal to n.
 */
int get_highest_power_of_2(int n){
    int res = 0;
    for (int i = n; i >= 1; i--) {
        // If i is a power of 2
        if ((i & (i - 1)) == 0) {
            res = i;
            break;
        }
    }
    return res;
}

/**
 * Converts a time domain signal into a frequency domain signal.
 * @param time_domain_signal The dignal to convert.
 * @param signal_size The size of the signal to convert.
 * @return An array containing the converted signal.
 */
void time_domain_to_frequency_domain(double time_domain_signal[], int signal_size, double frequency_domain_signal[], int new_signal_size){
    // Add real values to the output.
    for(int i = 0, j = 0; i < signal_size; i++, j += 2){
        frequency_domain_signal[j] = time_domain_signal[i];
    }

    // As per main() in 'test.c.'
    // Debugging note: read that frequency_domain_signal-1 could cause undefined behaviour.
    four1(frequency_domain_signal-1, new_signal_size/2, FFT_ISIGN);
}

/**
 * The main function of the program. Tests the convolve function with various input signals.
 * @param argc The number of command line arguments passed to the program. Should be greater than 2.
 * @param argv The command line arguments. Should be of the form <dry recording> <impulse response> <output file>
 * @return 0.
 */
int main(int argc, char* argv[]) {
    // Check for a valid number of command line arguments.
    if (argc != 4) {
        printf("Error: incorrect number of command line arguments.\n");
        return EXIT_FAILURE;
    }

    // Determine endianness of system.
    set_endianness();

    // Store the file names passed from the command line.
    std::string dry_recording_file = argv[1];
    std::string impulse_response_file = argv[2];
    std::string output_file = argv[3];

    // printf("dry_recording_file = %s\n", dry_recording_file.c_str());
    // printf("impulse_response_file = %s\n", impulse_response_file.c_str());
    // printf("output_file = %s\n", output_file.c_str());

    // Read input files.
    wav_file input_wav_file = read_wav_file(dry_recording_file);
    wav_file impulse_wav_file = read_wav_file(impulse_response_file);
    int output_file_size = input_wav_file.signal_size + impulse_wav_file.signal_size - 1;
    // print_wav_file(input_wav_file, "input_wav_file");
    // printf("\n");
    // print_wav_file(impulse_wav_file, "impulse_wav_file");

    // Start of FFT algorithm from class:
    // Convert time domain signals to frequency domain signals.
    double* x_time_domain = new double[input_wav_file.signal_size];
    double* h_time_domain = new double[impulse_wav_file.signal_size];
    prepare_signal_for_convolving(input_wav_file.signal, x_time_domain, input_wav_file.signal_size);
    prepare_signal_for_convolving(impulse_wav_file.signal, h_time_domain, impulse_wav_file.signal_size);
    int frequency_array_size = get_highest_power_of_2(std::max(input_wav_file.signal_size, impulse_wav_file.signal_size)) * 4; // Value must be large enough to prevent wrap around.
    double* x_frequency_domain = new double[frequency_array_size];
    double* h_frequency_domain = new double[frequency_array_size];
    for (int i = 0; i < frequency_array_size; i++) {
        x_frequency_domain[i] = 0.0;
        h_frequency_domain[i] = 0.0;
    }
    time_domain_to_frequency_domain(x_time_domain, input_wav_file.signal_size, x_frequency_domain, frequency_array_size);
    time_domain_to_frequency_domain(h_time_domain, impulse_wav_file.signal_size, h_frequency_domain, frequency_array_size);

    // Multiply x[k] by h[k] point by point.
    double* product = new double[frequency_array_size];
    for(int i = 0; i < frequency_array_size; i += 4){
        int complex_index_1 = i + 1;
        int complex_index_2 = i + 3;
        
        // First real part.
        product[i] = (x_frequency_domain[i] * h_frequency_domain[i]) - (x_frequency_domain[complex_index_1] * h_frequency_domain[complex_index_1]);

        // First complex part.
        product[complex_index_1] = (x_frequency_domain[complex_index_1] * h_frequency_domain[i]) + (x_frequency_domain[i] * h_frequency_domain[complex_index_1]);

        // Second real part.
        product[i + 2] = (x_frequency_domain[i + 2] * h_frequency_domain[i + 2]) - (x_frequency_domain[complex_index_2] * h_frequency_domain[complex_index_2]);

        // Second complex part.
        product[complex_index_2] = (x_frequency_domain[complex_index_2] * h_frequency_domain[i + 2]) + (x_frequency_domain[i + 2] * h_frequency_domain[complex_index_2]);
    }

    // Converting the result back to the time domain using the inverse FFT (IFFT).
    four1(product-1, frequency_array_size/2, INVERSE_FFT_ISIGN);
    double* convolved_signal = new double[output_file_size];
    for(int i = 0; i < output_file_size; i++){
        convolved_signal[i] = product[i * 2];
    }
    set_convolved_signal_magnitude(&input_wav_file, convolved_signal, output_file_size);

    // Prepare the output file and output file parameters.
    FILE* output_file_updated = fopen(output_file.c_str(), "wb");
    int input_file_bytes_per_sample = input_wav_file.bits_per_sample / BITS_PER_BYTE;
    short output_file_block_align = input_wav_file.num_channels * input_file_bytes_per_sample; 
    int output_file_subchunk_2_size = input_wav_file.num_channels * output_file_size * input_file_bytes_per_sample;
    int output_file_byte_rate = (int) input_wav_file.sample_rate * output_file_block_align;
    int output_file_chunk_size = output_file_subchunk_2_size + CHUNK_1_SIZE;

    // Write chunk_id.
    fputs(CHUNK1ID_DEFAULT, output_file_updated);

    // Write chunk_size.
    int size_of_int = sizeof(int);
    const size_t size_of_unsigned_char = sizeof(unsigned char);
    unsigned char formatted_int[size_of_int];
    memset(formatted_int, '\0', sizeof(formatted_int));
    format_int_endianness(output_file_chunk_size, formatted_int);
    fwrite(formatted_int, size_of_unsigned_char, size_of_int, output_file_updated);

    // Write format.
    fputs(FILE_TYPE_HEADER_DEFAULT, output_file_updated);

    // Write subchunk_1_id.
    fputs(SUBCHUNK_1_ID_DEFAULT, output_file_updated);

    // Write subchunk_1_size.
    memset(formatted_int, '\0', sizeof(formatted_int));
    format_int_endianness(SUBCHUNK_1_SIZE_NORMAL_VALUE, formatted_int);
    fwrite(formatted_int, size_of_unsigned_char, size_of_int, output_file_updated);

    // Write audio_format.
    int size_of_short = sizeof(short);
    unsigned char formatted_short[size_of_short];
    format_short_endianness(PCM_AUDIO_FORMAT_IDENTIFIER, formatted_short);
    fwrite(formatted_short, size_of_unsigned_char, size_of_short, output_file_updated);

    // Write num_channels.
    memset(formatted_short, '\0', sizeof(formatted_short));
    format_short_endianness(input_wav_file.num_channels, formatted_short);
    fwrite(formatted_short, size_of_unsigned_char, size_of_short, output_file_updated);

    // Write sample_rate.
    memset(formatted_int, '\0', sizeof(formatted_int));
    format_int_endianness(input_wav_file.sample_rate, formatted_int);
    fwrite(formatted_int, size_of_unsigned_char, size_of_int, output_file_updated);

    // Write byte_rate.
    memset(formatted_int, '\0', sizeof(formatted_int));
    format_int_endianness(input_wav_file.byte_rate, formatted_int);
    fwrite(formatted_int, size_of_unsigned_char, size_of_int, output_file_updated);

    // Write block_align.
    memset(formatted_short, '\0', sizeof(formatted_short));
    format_short_endianness(output_file_block_align, formatted_short);
    fwrite(formatted_short, size_of_unsigned_char, size_of_short, output_file_updated);

    // Write bits_per_sample.
    memset(formatted_short, '\0', sizeof(formatted_short));
    format_short_endianness(input_wav_file.bits_per_sample, formatted_short);
    fwrite(formatted_short, size_of_unsigned_char, size_of_short, output_file_updated);

    // Write subchunk_2_id.
    fputs(SUBCHUNK_2_ID_DEFAULT, output_file_updated);

    // Write subchunk_2_size.
    memset(formatted_int, '\0', sizeof(formatted_int));
    format_int_endianness(output_file_subchunk_2_size, formatted_int);
    fwrite(formatted_int, size_of_unsigned_char, size_of_int, output_file_updated);

    // Write convolved output signal into the output file.
    for(int i = 0; i < output_file_size; i++){
        memset(formatted_short, '\0', sizeof(formatted_short));
        format_short_endianness((short)convolved_signal[i], formatted_short);
        fwrite(formatted_short, size_of_unsigned_char, size_of_short, output_file_updated);
    }

    fclose(output_file_updated);

    return 0;
}