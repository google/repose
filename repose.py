#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 11:45:32 2020

@author: justinphillips
"""


# Process Holter V1.0 BETA - Utility for processing Holter signals
#
# Copyright 2020 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import datetime
import enum
import numpy as np
import pandas as pd
from pyedflib import highlevel
import os
from scipy import signal

 
class Pos(enum.Enum): 
  """ enum class for body position descriptive labels """
  SUPINE = 0
  PRONE = 1
  TILT = 2
  UPRIGHT = 3
  LEFT_SIDE = 4
  RIGHT_SIDE = 5


def get_data(pathname, filename):
  """reads data from files and returns dataframes for ecg and excel."""
  signals, signal_headers, header = highlevel.read_edf(os.path.join(pathname,
                                                                    filename))
  start_date = header['startdate']
  sample_rate = {header['label']:header['sample_rate']
                 for header in signal_headers}


  ecg_timestamp = [start_date +
                   datetime.timedelta(microseconds=1./sample_rate['ECG']*n*1e6)
                   for n in range(len(signals[0]))]  # 0 is index of ECG
  acc_timestamp = [start_date + datetime.timedelta
                   (microseconds=1./sample_rate['Accelerometer_X']*n*1e6)
                   for n in range(len(signals[1]))]  # 1 is index of Accel_X
  ecg = pd.DataFrame(signals[0], index=ecg_timestamp, columns=['ECG'])
  acc = pd.DataFrame(np.transpose(signals[1:4]), index=acc_timestamp,
            columns=[header['label'] for header in signal_headers[1:4]])
  acc['mag'] = np.sqrt(np.square(acc.Accelerometer_X)
                       + np.square(acc.Accelerometer_Y)
                       + np.square(acc.Accelerometer_Z))
  acc['total'] = acc['Accelerometer_X'] + acc['Accelerometer_Y'] + acc[
      'Accelerometer_Z']
  grav = np.median(acc.mag)
  acc.Accelerometer_X = acc.Accelerometer_X/grav
  acc.Accelerometer_Y = acc.Accelerometer_Y/grav
  acc.Accelerometer_Z = acc.Accelerometer_Z/grav
  acc.dropna()
  return ecg, acc, sample_rate


def butter_bandpass(lowcut, highcut, fs, order=2):
  """set params of bandpass filter for ECG signals."""
  nyq = 0.5 * fs
  low = lowcut / nyq
  high = highcut / nyq
  b, a = signal.butter(order, [low, high], btype='band')
  return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=2):
  """bandpass filter for ECG signals."""
  b, a = butter_bandpass(lowcut, highcut, fs, order=order)
  y = signal.lfilter(b, a, data)
  return y


def spectrum(sig, sample_rate):
  """returns amplitude spectrum of accel magnitude for RR."""
  mag = np.fft.rfft(sig)
  freq = np.fft.rfftfreq(len(sig), 1.0 / sample_rate)
  peakf = np.argmax(np.absolute(mag)) * freq[1]
  peakv = np.max(np.absolute(mag))
  return np.absolute(mag) / 1000., freq, peakf, peakv / 1000.


def calculate_hr(ecg, sample_rate, window_length=15, calc_every=5):
  """calculates heart rate from ECG."""
  ecg = ecg.resample('0.01S').mean()  # 100 Hz
  sample_rate['ECG'] = 100
  ecg_array = np.array(ecg.ECG)
  low, high = 3.0, 33.3   # filt. passband: remove baseline drift and HF noise
  _, _ = butter_bandpass(low, high, sample_rate['ECG'], order=2)
  ecg_filt = butter_bandpass_filter(ecg_array, low, high, sample_rate['ECG'],
                                    order=2)
  ecg_processed = -np.diff(ecg_filt)
  np.append(ecg_processed, 0)
  hr_win = window_length * sample_rate['ECG']
  hr = pd.DataFrame(columns=['HR'])
  for i in range(0, len(ecg_processed)-hr_win, sample_rate['ECG']*calc_every):
    ecg_sample = ecg_processed[i:i+hr_win]
    thresh = np.percentile(ecg_sample, 99)*0.75
    beats = []
    for j in range(len(ecg_sample)-1):
      if ecg_sample[j] < thresh < ecg_sample[j+1]:
        beats.append(j)
    hr_mean = 60*1/(np.mean(np.diff(beats)/sample_rate['ECG']))
    new_row = pd.DataFrame([[hr_mean]], columns=['HR'],
                           index=[ecg.index[i+hr_win]])
    hr = pd.concat([hr, pd.DataFrame(new_row)], ignore_index=False)
    hr = hr.resample('S').pad()
  return hr


def calculate_rr(acc, sample_rate, window_length=25, calc_every=5):
  """calculates respiratory rate from accelerometer sum."""
  acc_tot = np.array(acc.total)  # x + y + z
  rr_win = window_length * sample_rate['Accelerometer_X']
  rr = pd.DataFrame(columns=['RR'])
  for i in range(0, len(acc_tot)-rr_win, sample_rate['ECG']*calc_every):
    windowed_acc = np.subtract(acc_tot[i:i+rr_win],
                               np.median(acc_tot[i:i+rr_win]))
    _, _, rf, _ = spectrum(windowed_acc, sample_rate['Accelerometer_X'])
    new_row = pd.DataFrame([[rf*60]], columns=['RR'],
                           index=[acc.index[i+rr_win]])
    rr = pd.concat([rr, pd.DataFrame(new_row)], ignore_index=False)
    rr = rr.resample('S').pad()
  return rr


def body_pos_and_angles(acc):
  """calculates body position (as index# and description)."""

  time, pos, pos_name, thetas, phis = [], [], [], [], []
  for i, _ in acc.iterrows():
    acc_x = np.clip(acc.Accelerometer_X[i], -1., 1.)
    acc_y = np.clip(acc.Accelerometer_Y[i], -1., 1.)
    acc_z = np.clip(acc.Accelerometer_Z[i], -1., 1.)
    time.append(i)
    theta = np.arcsin(acc_x)*180/np.pi
    phi = np.sign(np.arcsin(-acc_y))*np.arccos(-acc_z)*180/np.pi
    # decision trees for body position (affects angles calculations)
    if abs(acc_y) > 0.5:  # y >+- 30˚
      if acc_y < 0:
        position = Pos.LEFT_SIDE
      else:
        position = Pos.RIGHT_SIDE
      thetas.append(float('nan'))
    else:
      if acc_z > acc_x and acc_z > acc_y:
        position = Pos.PRONE
      elif acc_z < -0.966:   # z < 15˚
        position = Pos.SUPINE
      elif acc_x > 0.966:    # x > 75˚
        position = Pos.UPRIGHT
      else:
        position = Pos.TILT
      thetas.append(theta)

    if acc_x > 0.71:
      phis.append(float('nan'))
    else:
      phis.append(phi)

    pos.append(position.value)
    pos_name.append(position.name)
  return pos, pos_name, thetas, phis


def main():
  pathname = '//Users/<user-name>/Documents/EDF_Files/20200513'  # <-- * 
  # (* folder name (containing EDFs) goes here)
  for filename in os.listdir(pathname):
    if '.EDF' not in filename and '.edf' not in filename:
      continue
    print('Reading:', os.path.join(pathname, filename))
    ecg, acc, sample_rate = get_data(pathname, filename)

    print('Recording start:', ecg.index[0].ctime())
    print('Recording end:  ', ecg.index[-1].ctime())

    hr = calculate_hr(ecg, sample_rate, window_length=15, calc_every=5) # <- *
    rr = calculate_rr(acc, sample_rate, window_length=25, calc_every=5) # <- *
    # (* averaging window and 'calculate every x seconds' for HR and RR)
    
    acc = acc.resample('S').mean()   # <-- * downsample outputs to 1 S/s
    # (* note that the accelerometer data is initially 
    #    downsampled to this freq then all other outputs are synced to this.)
    
    pos, pos_name, theta, phi = body_pos_and_angles(acc)

    out = pd.DataFrame(index=acc.index)
    out = pd.concat([out, rr], ignore_index=False, axis=1)
    out = pd.concat([out, hr], ignore_index=False, axis=1)
    out['Pos_Index'] = pos
    out['Position'] = pos_name
    out['Tilt'] = theta
    out['Rotation'] = phi

    excel_filename = filename.strip('.EDF')+'.xlsx'
    print('Writing:', os.path.join(pathname, excel_filename))

          
if __name__ == '__main__':
  main()
    
    
   
