import matplotlib.pyplot as plt
import numpy as np
import excel_reader as exr


def compareForceModelData(fileControl, filePathPerturbation):
       # output rMag array
       # Initialize the reader
       readerControl = exr.ExcelReader(fileControl)
       readerPerturb= exr.ExcelReader(filePathPerturbation)

       # Read first sheet
       df_control = readerControl.read_sheet()
       df_perturb = readerPerturb.read_sheet()

       tVec_control = df_control["Time(s)"]
       r_control = df_control[["X(km)","Y(km)","Z(km)"]]

       tVec_perturb = df_perturb["Time(s)"]
       r_perturb = df_perturb[["X(km)","Y(km)","Z(km)"]]

       # Option 1: Convert to numpy array first (more efficient)
       r_control_array = r_control.to_numpy()
       r_perturb_array = r_perturb.to_numpy()
       rDiffMag = np.full(len(r_control), np.nan)
       for i in range(len(tVec_control)):
              if tVec_control[i] != tVec_perturb[i]:
                     raise ValueError("Whoops")
              rDiff = r_control_array[i] - r_perturb_array[i]
              rDiffMag[i] = np.linalg.norm(rDiff)

       return tVec_control, rDiffMag


dataTwoBody = "/Users/douglasgarza/CLionProjects/atlas/test_data/force_model_magnitude/two_body_pos_data.csv"
data2x0 = "/Users/douglasgarza/CLionProjects/atlas/test_data/force_model_magnitude/2x0_pos_data.csv"
data3x0 = "/Users/douglasgarza/CLionProjects/atlas/test_data/force_model_magnitude/3x0_pos_data.csv"
dataDrag = "/Users/douglasgarza/CLionProjects/atlas/test_data/force_model_magnitude/drag_pos_data.csv"
dataSRP = "/Users/douglasgarza/CLionProjects/atlas/test_data/force_model_magnitude/srp_pos_data.csv"
data3Body = "/Users/douglasgarza/CLionProjects/atlas/test_data/force_model_magnitude/thirdBody_pos_data.csv"


tVec_twoBody, rDiffMag2x0 = compareForceModelData(dataTwoBody, data2x0)
tVec_2x0, rDiffMag3x0 = compareForceModelData(data2x0, data3x0)
tVec_twoBody, rDiffMagDrag = compareForceModelData(dataTwoBody, dataDrag)
tVec_twoBody, rDiffMagSRP = compareForceModelData(dataTwoBody, dataSRP)
tVec_twoBody, rDiffMag3Body = compareForceModelData(dataTwoBody, data3Body)

xtick_positions = np.array([0, 1440, 2880, 4320, 5760])
xtick_labels = ['0', '1440', '2880', '4320', '5760']
# ytick_positions = np.array([0.1, 1.0, 10, 100, 1000, 1000, 10000, 100000, 1000000])
# ytick_labels = ['0.1', '1.0', '10', '100', '1000', '1000', '10000', '100000', '1000000']

fig, ax = plt.subplots()
ax.plot(tVec_twoBody / 60, rDiffMag2x0 * 1000)
ax.plot(tVec_twoBody / 60, rDiffMagDrag * 1000)
ax.plot(tVec_twoBody / 60, rDiffMagSRP * 1000)
ax.plot(tVec_twoBody / 60, rDiffMag3Body * 1000)


ax.set(xlabel='time (min)', ylabel='Difference (m)',
       title='Force Model Comparisons')
plt.xticks(xtick_positions, xtick_labels)
# plt.yticks(ytick_positions, ytick_labels)
ax.grid()     
ax.set_yscale('log')  
plt.show()






