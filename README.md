# Name

Geometry-based kinetic simulation of DNA-gold nanoparticle motor
https://github.com/TakanoriHarashima/DNA_NP_motor_simulation/blob/main/DNAmotor_simulation_Gillespie_v03.py

# Manuscript

Takanori Harashima, Akihiro Otomo, and Ryota Iino, Rational engineering of DNA-nanoparticle motor with high speed and processivity comparable to motor proteins (under revision)

# Features

Details of the algorithm and description of the simulation parameters are referred in the manuscript (https://www.biorxiv.org/content/10.1101/2024.05.23.595615v2). First, a two-dimensional pixel matrix was defined to model the RNAs on the glass surface. Pixel size was normalized by the DNA density on the AuNP. RNAs were randomly distributed with a ratio of the RNA density to the DNA density. Reaction at a single RNA site is assumed to contain three sequential elementary steps: DNA/RNA hybridization, RNase H binding to DNA/RNA duplex, and RNA hydrolysis. The reaction proceeds within the accessible area of DNA with the radius of 28.0 nm. Rate constants of the three elementary steps, konDNA/RNA, k E, and kcatE are the simulation parameters. Simulation steps proceeded in units of reaction events. ***The reaction site and time step is determined based on Gillespie algorithm.*** Constraint of the motor position was considered by introducing mobile region of each DNA/RNA hybrid (25.2 nm). The motor position was determined randomly within the region where all mobile regions overlap. Unless all mobile regions overlapped, the motor position was fixed. 

# Environment Tested
  - Windows 10 Pro  
  - Anaconda 3  
  - python 3.9  
  - Spyder IDE 5.4.1  
  - ffmpeg 4.2.3  

(Non-standard hardware is not required.)  

# Requirement
Python packages
| Package  | Version |
| ------------- | ------------- |
| pandas  | 1.4.4  |
| numpy  | 1.21.5  |
| scipy  | 1.9.1  |
| opencv-python  | 4.7.0.68  |
| matplotlib  | 3.5.2  |
| psutil  | 5.9.0  |
| glob2  | 0.7  |

# Installation Guide
Whole installation procedure typically takes ~1 hour on a normal desktop computer.  
**1. Install Anaconda:**  
  - https://www.anaconda.com/download

**2. Open Spyder:**  
  - To run the bundled version of Spyder after installing it with Anaconda, the recommended method on Windows is to launch it via the Start menu shortcut. On other platforms, open Anaconda Navigator, scroll to Spyder under Home and click Launch.
  - https://docs.spyder-ide.org/current/installation.html

**3. Install required python modules:**  
  - Download requirements.txt from Github.  
  - Open a console and install required python modules.  
    ```pip install -r requirements.txt```
    
**4. Run the simulation:**  
  - Download 	***DNAmotor_simulation_Gillespie_v03.py*** from Github.
  - Spyder -> File -> Open -> Select DNAmotor_simulation_Gillespie_v03.py
  - Set parameters.
  - Run File (push F5 key) and simulation starts.

# Demo
**1. Setting:**  
See the manuscript and values of parameters are described in Supplementary Table 3.  
Just for a demo, here is the recommended parameters  
```
# Basic parameters
globalfol = r'Directory\to\perform\DNAmotor\simulation'
date = 'DATE'
RNA_size = 2000                        # Full range of RNA substrate (px)  (default:2000)
tmax_simu = 100000                     # Uplimit for simulation (sec)  (default:100000)
N_simu = 3                             # Nunber of trajectory generate (default:5)
frame_per_event = 1000                 # Span frame to check the progress of simulation (default:1000)
foli=0                                 # ID of the condition of kinetic parameters (default:0)
# Kinetic parameters
khyb_list = np.array([0.2])                  # DNA/RNA Hybridization rate [s-1]　
konE_list = np.array([1.0]) *10**6           # RNase H binding rate [M-1 s-1]
kcat_list = np.array([3])                    # RNA hydrolysis rate [s-1] 
RNaseH_list = np.array([36])                 # RNase H condition [nM]
#######################
SIMULATION = True # Run the new simulation
SUMMERIZE = True # Output rough summary of the simulation
```
**2. Run the simulation:**  
Open DNAmotor_simu_v5.03.py with Spyder.  
Push F5 key to run the simulation.  

**3. Expected output:**  
The simulation make folders as follow.  
<pre>
001_khyb=0.30_kcatE=4.0_konE=1.0x106.
└─DATE_{'N_simu'= 3, 'tmax'= 100000, 'RNaseH'= 36, 'frame_per_event'= 1000}
    └─progress
        ├─000
        ├─001
        └─002
</pre>
You can see the progress of the simulation by looking the output figure files which will be generated in individual directories (i.e. '000', '001', '002').  
Expected run time for demo on a normal desktop computer is 2 min/trajectory.  

**4. Make resulting movies:**  
Download DNAmotor_simu_movie_maker_per_time_v05.py from Github and open with Spyder.
Set the parameter. Here is an example.
```
workfol = r'Directory\to\perform\DNAmotor\simulation\001_khyb=0.30_kcatE=4.0_konE=1.0x106'
path_list = [workfol + os.sep + r"DATE_{'N_simu'= 3, 'tmax'= 100000, 'RNaseH'= 36, 'frame_per_event'= 1000}"]
fps = 20                         # Frame per second, adjust to experimental conditions
dt_timecourse = 500              # Range to show time-course [sec]
itrace= 0                        # Trajectory ID to make the movie 
###############################################################################
###############################################################################
MAKE_FIGURES = True # If True, the program generate every snapshots during simulation as png data.  
MAKE_MOVIE = True   # If True, the program concatenate the images by ffmpeg.  
```  
This program make the image of snapshot and concatenate them into the mp4 movie (see Supplementary Movie 4,5,6,7).

# Reproduction instructions  
Here is the parameter set for our simulation shown in Figure ***3***.  
```
# Basic parameters
globalfol = r'Directory\to\perform\DNAmotor\simulation'
date = 'DATE'
RNA_size = 3000                        # Full range of RNA substrate (px)  (default:2000)
tmax_simu = 100000                     # Uplimit for simulation (sec)  (default:100000)
N_simu = 50                             # Nunber of trajectory generate 
frame_per_event = 1000                 # Span frame to check the progress of simulation (default:1000)
foli=0                                 # ID of the condition of kinetic parameters (default:0)
# Kinetic parameters
khyb_list = np.array([0.2])                  # DNA/RNA Hybridization rate [s-1]　
konE_list = np.array([1.0]) *10**6           # RNase H binding rate [M-1 s-1]
kcat_list = np.array([3])                    # RNA hydrolysis rate [s-1] 
RNaseH_list = np.array([36,144,360,720,1440,3600,7200])                 # RNase H condition [nM]
#######################
SIMULATION = True # Run the new simulation
SUMMERIZE = True # Output rough summary of the simulation
```
  
Here is the parameter set for the simulation-based fitting shown in Supplementary Fig. 11-16.
```
# Basic parameters
globalfol = r'Directory\to\perform\DNAmotor\simulation'
date = 'DATE'
RNA_size = 3000                        # Full range of RNA substrate (px)  (default:2000)
tmax_simu = 100000                     # Uplimit for simulation (sec)  (default:100000)
N_simu = 20                             # Nunber of trajectory generate 
frame_per_event = 1000                 # Span frame to check the progress of simulation (default:1000)
foli=0                                 # ID of the condition of kinetic parameters (default:0)
# Kinetic parameters
khyb_list = np.array([0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 5.0])                  # DNA/RNA Hybridization rate [s-1]
konE_list = np.array([0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 5.0]) *10**6           # RNase H binding rate [M-1 s-1]
kcat_list = np.array([0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0])                    # RNA hydrolysis rate [s-1] 
RNaseH_list = np.array([36,144,360,720,1440,3600,7200])                 # RNase H condition [nM]
#######################
SIMULATION = True # Run the new simulation
SUMMERIZE = True # Output rough summary of the simulation
```


# Author

* Dr. Takanori Harashima
* Institute for Molecular Science, National Institutes of Natural Sciences, Okazaki, Aichi 444-8787, Japan
* harashima@ims.ac.jp

# License

Our simulation code used in this study is provided under the MIT License.

