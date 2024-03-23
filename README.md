# Convolution Reverb

This program applies reverb to an audio file by combining it with a given impulse response.
Audio signal convolution is done using the overlap-add Fast Fourier Transform algorithm. 

## Usage

To compile the program, use g++ or any other C++ compiler:

```bash
g++ convolve.cpp -o convolve
```

Once compiled, you can run the program with the following command:

```bash
./convolve <dry_recording.wav> <impulse_response.wav> <output_file.wav>
```
- **<dry_recording.wav>:** Specifies the name of a dry recording WAV file.
- **<impulse_response.wav>:** Specifies the name of an impulse response WAV file.
- **<output_file.wav>:** Specifies the name of the output WAV file.
