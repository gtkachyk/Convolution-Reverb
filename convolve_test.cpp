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
#define BITS_PER_BYTE 8
#define STRING_FIELD_SIZE 4
#define SUBCHUNK_1_SIZE_NORMAL_VALUE 16
#define BITS_PER_SAMPLE_DEFAULT 16

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
 * Compares two char arrays to determine if they are equal.
 * @param wav_file_1 The first char array of the equality.
 * @param wav_file_2 The second char array of the equality.
 * @param size The size of the arrays.
 * @return True if array_1 and array_2 have all the same element values, false otherwise.
 */
bool char_arrays_are_arrays_equal(const char array_1[], const char array_2[], int size) {
    for (int i = 0; i < size; i++) {
        // printf("array_1[%d] = %d, array_2[%d] = %d\n", i, array_1[i], i, array_2[i]);
        if (array_1[i] != array_2[i]) {
            printf("array_1[%d] = %d, array_2[%d] = %d\n", i, array_1[i], i, array_2[i]);
            return false;
        }
    }
    return true;
}

/**
 * Compares two short arrays to determine if they are equal.
 * @param wav_file_1 The first short array of the equality.
 * @param wav_file_2 The second short array of the equality.
 * @param size The size of the arrays.
 * @return True if array_1 and array_2 have all the same element values, false otherwise.
 */
bool short_arrays_are_arrays_equal(const short array_1[], const short array_2[], int size) {
    for (int i = 0; i < size; i++) {
        // printf("array_1[%d] = %d, array_2[%d] = %d\n", i, array_1[i], i, array_2[i]);
        if (array_1[i] != array_2[i]) {
            printf("array_1[%d] = %d, array_2[%d] = %d\n", i, array_1[i], i, array_2[i]);
            return false;
        }
    }
    return true;
}

/**
 * Compares the parameters of two wav_file structs to determine if they are equal.
 * @param wav_file_1 The first wav_file struct of the equality.
 * @param wav_file_2 The second wav_file struct of the equality.
 * @return True if wav_file_1 and wav_file_2 have all the same parameter values, false otherwise.
 */
bool are_wav_files_equal(wav_file wav_file_1, wav_file wav_file_2){
    // RIFF chunk descriptor.
    if(!char_arrays_are_arrays_equal(wav_file_1.chunk_id, wav_file_2.chunk_id, sizeof(wav_file_1.chunk_id))){
        printf("The wav file have different chunk_id.\n");
        return false;
    }
    if(wav_file_1.chunk_size != wav_file_2.chunk_size){
        printf("The wav file have different chunk_size.\n");
        return false;
    }
    if(!char_arrays_are_arrays_equal(wav_file_1.format, wav_file_2.format, sizeof(wav_file_1.format))){
        printf("The wav file have different format.\n");
        return false;
    }

    // fmt sub-chunk.
    if(!char_arrays_are_arrays_equal(wav_file_1.subchunk_1_id, wav_file_2.subchunk_1_id, sizeof(wav_file_1.subchunk_1_id))){
        printf("The wav file have different subchunk_1_id.\n");
        return false;
    }
    if(wav_file_1.subchunk_1_size != wav_file_2.subchunk_1_size){
        printf("The wav file have different subchunk_1_size.\n");
        return false;
    }
    if(wav_file_1.audio_format != wav_file_2.audio_format){
        printf("The wav file have different audio_format.\n");
        return false;
    }
    if(wav_file_1.num_channels != wav_file_2.num_channels){
        printf("The wav file have different num_channels.\n");
        return false;
    }
    if(wav_file_1.sample_rate != wav_file_2.sample_rate){
        printf("The wav file have different sample_rate.\n");
        return false;
    }
    if(wav_file_1.byte_rate != wav_file_2.byte_rate){
        printf("The wav file have different byte_rate.\n");
        return false;
    }
    if(wav_file_1.block_align != wav_file_2.block_align){
        printf("The wav file have different block_align.\n");
        return false;
    }
    if(wav_file_1.bits_per_sample != wav_file_2.bits_per_sample){
        printf("The wav file have different bits_per_sample.\n");
        return false;
    }

    // data sub-chunk.
    if(!char_arrays_are_arrays_equal(wav_file_1.subchunk_2_id, wav_file_2.subchunk_2_id, sizeof(wav_file_1.subchunk_2_id))){
        printf("The wav file have different subchunk_2_id.\n");
        return false;
    }
    if(wav_file_1.subchunk_2_size != wav_file_2.subchunk_2_size){
        printf("The wav file have different subchunk_2_size.\n");
        return false;
    }
    if(wav_file_1.subchunk_2_size == wav_file_2.subchunk_2_size){
        if(!char_arrays_are_arrays_equal(wav_file_1.data, wav_file_2.data, wav_file_1.subchunk_2_size)){
            printf("The wav file have different data arrays.\n");
            return false;
        }  
    }
    else{
        printf("Cannot compare arrays of different lengths.\n");
        return false;
    }
    if(wav_file_1.signal_size == wav_file_2.signal_size){
        if(!short_arrays_are_arrays_equal(wav_file_1.signal, wav_file_2.signal, wav_file_1.signal_size)){
            printf("The wav file have different signal arrays.\n");
            return false;
        }   
    }
    else{
        printf("Cannot compare arrays of different lengths.\n");
        return false;
    }
    if(wav_file_1.signal_size != wav_file_2.signal_size){
        printf("The wav file have different signal_size\n");
        return false;   
    }
    return true;
}

/**
 * The main function of the program. Tests the convolve function with various input signals.
 * @param argc The number of command line arguments passed to the program. Should be greater than 2.
 * @param argv The command line arguments. Should be of the form <dry recording> <impulse response> <output file>
 * @return 0.
 */
int main(int argc, char* argv[]) {
    // Check for a valid number of command line arguments.
    if (argc != 3) {
        printf("Error: incorrect number of command line arguments.\n");
        return EXIT_FAILURE;
    }

    // Store the file names passed from the command line.
    std::string file_1 = argv[1];
    std::string file_2 = argv[2];

    printf("file_1 = %s\n", file_1.c_str());
    printf("file_2 = %s\n", file_2.c_str());

    // Read input files.
    wav_file file_1_wav_file = read_wav_file(file_1);
    wav_file file_2_wav_file = read_wav_file(file_2);
    print_wav_file(file_1_wav_file, "file_1_wav_file");
    printf("\n");
    print_wav_file(file_2_wav_file, "file_2_wav_file");
    
    if(are_wav_files_equal(file_1_wav_file, file_2_wav_file)){
        printf("The wav files are equal.\n");
    }
    else{
        printf("The wav files are not equal.\n");
    }


    return 0;
}