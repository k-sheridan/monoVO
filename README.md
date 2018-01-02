## mono VO
Developed by Kevin Sheridan, Purdue University.

While this software does work, it is lacking occlusion detection and a bugless EKF implementation.

I extracted this code from my INVIO project because I decided to switch my method for visual odometry. This visual odometry algorithm is not robust enough for my liking due to how it estimates the depth of features and deals with outliers. 


## Performance

![Small Scale Results](/images/monoVO.png)

![Brick Wall Results](/images/monoVO2.png)

The speed of the program per frame on a laptop (2015 Macbook Pro in my case) 
### runtimes per frame:
Feature Extraction: 0.1-0.3ms 

KLT Feature Tracking for 300 features: 0.8-3ms 

Motion Estimation with 100-300 features: 0.03 - 0.2ms 

Depth Measurement and Update for 200 features: 0.5ms

###Total: 1.43 - 5.43ms 

Runs using 1 thread.

## Quick ROS installation guide

(developed using ROS Kinetic)

1. >cd ~/catkin_ws/src
2. >git clone https://github.com/pauvsi/monoVO
3. >sudo apt-get install ros-{your distribution}-sophus
4. >cd ..
5. >catkin_make

this should compile the entire package. If there is an issue please create an issue on this repo!

## Usage

You must provide a rectified mono image with a corresponding camera info on topics /camera/image_rect and /camera/camera_info respectively. This is temporary.

## Development Status

Finished for the foreseeable future.

