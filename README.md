# Interpolated Cross-Correlation Function (ICCF) Analysis

Image from python code:
<img src="https://github.com/PuDu-Astro/images_for_doc/blob/master/ICCF_python.png" width="720" height="438">

Image from c code:
<img src="https://github.com/PuDu-Astro/images_for_doc/blob/master/ICCF_c.png" width="720" height="438">

Interpolated Cross-Correlation Function (ICCF) analysis for light curves, with Monte Carlo error estimation.

**Author**: Pu Du  
**Affiliation**: Institute of High Energy Physics (IHEP)  
**Email**: dupu@ihep.ac.cn  

## 📖 Reference
This implementation is based on:
- Gaskell & Sparke (1986), ApJ, 305, 175 (core algorithm)
- Peterson et al. (1998), PASP, 110, 660 (FR/RSS error estimation method)

## 🛠 Installation & Dependencies

### System Requirements
- **C dependencies**:
  - [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
  - [PLplot](http://plplot.sourceforge.net/)

- **Python dependencies**:
  - numpy (`pip install numpy`)
  - matplotlib (`pip install matplotlib`)

## 🚀 Quick Start

### C Version
1. Compile the C program:
   ```bash
   gcc -O2 piccf_mc.c -o piccf_mc -lgsl -lplplot -lm

2. Run the executable:
   ```bash
   ./piccf_mc
   ```
   The program will read parameters from the param file.

### Python Version
1. First compile the C code as a shared library:
   ```bash
   gcc -O2 -fPIC -shared piccf_mc.c -o libpiccf_mc.so -lgsl -lplplot -lm
   ```

2. Update `LIB_PATH` in `piccf_mc.py` to point to your `.so` library
   
4. Run either:
   ```bash
   # Method 1: Direct parameter input
   python3 piccf_mc.py con_name line_name tbeg tend nt nbin num_mc
   
   # Method 2: Parameter file
   python3 piccf_mc.py param
   ```

## ⚙️ Configuration Parameters
|Parameter	|Example Value	|Description  |
|:----------|:--------------|:------------|
|continuum	|con.txt	      |Continuum light curve data file|
|line	      |line.txt	      |Emission-line light curve data file|
|tbeg	      |-10.0	        |Lower time limit for CCF analysis (days)|
|tend	      |120.0	        |Upper time limit for CCF analysis (days)|
|nt	        |1301	          |Number of time bins for CCF calculation|
|nbin	      |131	          |Number of bins for FR/RSS result display|
|num_mc	    |1000	          |Number of Monte Carlo simulations for FR/RSS|

## 📊 Output
The program generates:
- ICCF curve with peak lag identification
- Cross-correlation centroid distribution (CCCD) and cross-correlation peak distribution (CCPD) from FR/RSS
- Statistical summary of time lag measurements

### Output files:
|File	|Description  |
|:----------|:--------------|
|piccf_out.txt	|ICCF	      |
|peak_mc.txt	  |peak values from FR/RSS Monte Carlo |
|cent_mc.txt	  |centroid values from FR/RSS Monte Carlo	        |
|peak_dist.txt	|CCPD from FR/RSS	        |
|cent_disk.txt  |CCCD from FR/RSS	          |
|pypiccf_ccf_out.txt	  |ICCF (python version)	          |
|pypiccf_mc_out.txt	    |CCCD and CCPD (python version)	          |
