# example_mdanalysis

# Prerequisites & Installation

1. Virtual python environments via Mamba/Anaconda/Miniconda
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```

2. Create virtual environment
```
mamba create -n [env_name] -c conda-forge mdanalysis jupyter openmm=8.3 meson=1.1.1
mamba activate [env_name]
```

3. Compile libraries via f2py in Numpy
```
cd lib
python -m numpy.f2py -c -m mylib subroutines.f90
```

4. Usage
```
mamba activate [env_name]
# For Anaconada/Miniconda
# conda activate [env_name]

# Open ipynb files
```

5. How to use jupyter lab
```
ssh -L 8888:localhost:8888 [id]@[servername]
mamba activate [env_name]
cd example_mdanalysis
nohup jupyter lab --port [port number] &
tail nohup.out # and copy jupyter server link
# And open that link in your own local browser!
```
