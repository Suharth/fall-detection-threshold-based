# fall-detection-threshold-based
MATLAB implementation of the threshold based fall detection algorithm proposed in the paper:  
A.K. Bourke, G.M. Lyons,
A threshold-based fall-detection algorithm using a bi-axial gyroscope sensor,
Medical Engineering & Physics,
Volume 30, Issue 1,
2008,
Pages 84-90,
ISSN 1350-4533,
https://doi.org/10.1016/j.medengphy.2006.12.001.
(https://www.sciencedirect.com/science/article/pii/S1350453306002657)
Abstract: A threshold-based algorithm, to distinguish between Activities of Daily Living (ADL) and falls is described. A gyroscope based fall-detection sensor array is used. Using simulated-falls performed by young volunteers under supervised conditions onto crash mats and ADL performed by elderly subjects, the ability to discriminate between falls and ADL was achieved using a bi-axial gyroscope sensor mounted on the trunk, measuring pitch and roll angular velocities, and a threshold-based algorithm. Data analysis was performed using Matlab® to determine the angular accelerations, angular velocities and changes in trunk angle recorded, during eight different fall and ADL types. Three thresholds were identified so that a fall could be distinguished from an ADL: if the resultant angular velocity is greater than 3.1rads/s (Fall Threshold 1), the resultant angular acceleration is greater than 0.05rads/s2 (Fall Threshold 2), and the resultant change in trunk-angle is greater than 0.59rad (Fall Threshold 3), a fall is detected. Results show that falls can be distinguished from ADL with 100% accuracy, for a total data set of 480 movements.
Keywords: Falls in the elderly; Fall detection; Gyroscope; Activities of Daily Living; Threshold

Modify the initial part of the code to use any dataset that has roll and pitch angular velocity values. The dataset I used can be found in the following Github repository:
https://github.com/nhoyh/HR_IMU_falldetection_dataset
