clc
files = dir('/MATLAB Drive/ASEN 2002 Lab/Lab 2/Data collection/VTtoPT');
long_name = strcat(files(5).folder,'/',files(3).name);
data1 = load(long_name);

files2 = dir('/MATLAB Drive/ASEN 2002 Lab/Lab 2/Data collection/PPtoPT');
long_name2 = strcat(files2(6).folder,'/',files(7).name);
data2 = load(long_name2);