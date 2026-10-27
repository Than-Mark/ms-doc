import matplotlib.colors as mcolors
import matplotlib.animation as ma
import matplotlib.pyplot as plt
from typing import List, Tuple
from itertools import product
import numpy.typing as npt
import pandas as pd
import numpy as np
import numba as nb
import imageio
import sys
import os
import shutil

randomSeed = 10

if "ipykernel_launcher.py" in sys.argv[0]:
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm

colors = ["#5657A4", "#95D3A2", "#FFFFBF", "#F79051", "#A30644"]
cmap = mcolors.LinearSegmentedColormap.from_list("my_colormap", colors)
new_cmap = mcolors.LinearSegmentedColormap.from_list(
    "new", plt.cm.hsv(np.linspace(0, 1, 256)) * 0.85, N=256
)

@nb.njit
def colors_idx(phaseTheta):
    return np.floor(256 - phaseTheta / (2 * np.pi) * 256).astype(np.int32)


plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
if os.path.exists("/opt/conda/bin/ffmpeg"):
    plt.rcParams['animation.ffmpeg_path'] = "/opt/conda/bin/ffmpeg"
else:
    plt.rcParams['animation.ffmpeg_path'] = "D:/Programs/ffmpeg/bin/ffmpeg.exe"

sys.path.append("..")
from swarmalatorlib.template import Swarmalators2D


class StateAnalysis:
    def __init__(self, model):
        if model is None:
            return
        self.model = model
        
        targetPath = f"{self.model.savePath}/{self.model}.h5"
        
        totalPositionX = pd.read_hdf(targetPath, key="positionX")
        TNum = totalPositionX.shape[0] // self.model.agentsNum
        self.TNum = TNum

        self.totalPositionX = totalPositionX.values.reshape(TNum, self.model.agentsNum, 2)

    def get_state(self, index: int = -1):
        positionX = self.totalPositionX[index]

        return positionX

    def plot(self, index: int = -1, ax: plt.Axes = None, particleStyle: str = 'scatter'):
    
        assert particleStyle in ['scatter']

        if ax is None:
            _, ax = plt.subplots(figsize=(5, 5))

        positionX = self.get_state(index)

        if particleStyle == 'scatter':
            ax.scatter(positionX[:, 0], positionX[:, 1], s=3, color="black", linewidths=0)
        # elif particleStyle == 'quiver':
        #     norm = np.linalg.norm(self.temp["dotR"], axis=1, keepdims=True)
        #     norm[norm == 0] = 1
        #     unitDirection = self.temp["dotR"] / norm
            
        #     ax.quiver(
        #         self.positionX[:, 0], self.positionX[:, 1],
        #         unitDirection[:, 0], unitDirection[:, 1],
        #         angles="xy", scale_units="xy"
        #     )
        ax.set_xlim(0, self.model.boundaryLength)
        ax.set_ylim(0, self.model.boundaryLength)


class QuorumSensing(Swarmalators2D):
    stateAnalysisClass = StateAnalysis
    def __init__(self, distanceDc: float, 
                 thresLp: float = 1.0, boundaryLength: float = 10, 
                 diffusionDmax: float = 0.1, diffusionDmin: float = 0.01,
                 sensingStyle: str = "inverse_distance_sum",
                 agentsNum: int = 1000, dt: float = 0.01,
                 tqdm: bool = False, savePath: str = None, shotsnaps: int = 10,
                 randomSeed: int = 10, overWrite: bool = False) -> None:

        assert sensingStyle in ["neighbor_count", "inverse_distance_sum"]
        
        self.distanceDc = distanceDc
        self.thresLp = thresLp
        self.boundaryLength = boundaryLength
        self.diffusionDmax = diffusionDmax
        self.diffusionDmin = diffusionDmin

        self.agentsNum = agentsNum
        self.dt = dt
        self.tqdm = tqdm
        self.savePath = savePath
        self.shotsnaps = shotsnaps
        self.randomSeed = randomSeed
        self.overWrite = overWrite
        
        np.random.seed(randomSeed)
        self.positionX = np.random.random((agentsNum, 2)) * boundaryLength

        self.halfBoundaryLength = boundaryLength / 2
        self.rho0 = agentsNum / (boundaryLength ** 2)
        self.sensingArea = np.pi * distanceDc ** 2
        if sensingStyle == "neighbor_count":
            self.rhoThres = thresLp * self.rho0
        elif sensingStyle == "inverse_distance_sum":
            self.Pth = thresLp * self.rho0
        self.counts = 0
        self.sensingStyle = sensingStyle
        self.dotRParams = (
            self.boundaryLength,
            self.halfBoundaryLength,
            self.distanceDc,
        )

        self.temp = dict()
        self.temp["dotR"] = self.dotR

    @staticmethod
    @nb.njit
    def _neighbor_count_nI(positionX: np.ndarray, params: Tuple[float]) -> np.ndarray:
        agentsNum = positionX.shape[0]
        boundaryLength, halfBoundaryLength, distanceD0 = params

        nI = np.zeros(agentsNum)
        for i in range(agentsNum):
            xDiff = np.abs(positionX[:, 0] - positionX[i, 0])
            yDiff = np.abs(positionX[:, 1] - positionX[i, 1])
            neighborIdxs = np.where(
                (xDiff < distanceD0) | (boundaryLength - xDiff < distanceD0) & 
                (yDiff < distanceD0) | (boundaryLength - yDiff < distanceD0)
            )[0]
            if neighborIdxs.size == 0:
                continue

            subX = positionX[i] - positionX[neighborIdxs]
            deltaX = positionX[i] - (
                positionX[neighborIdxs] * (-halfBoundaryLength <= subX) * (subX <= halfBoundaryLength) + 
                (positionX[neighborIdxs] - boundaryLength) * (subX < -halfBoundaryLength) + 
                (positionX[neighborIdxs] + boundaryLength) * (subX > halfBoundaryLength)
            )
            distance = np.sqrt(np.sum(deltaX**2, axis=1))
            A = np.where(distance <= distanceD0)[0]

            nI[i] = len(A)

        return nI

    @staticmethod
    @nb.njit
    def _Pi(positionX: np.ndarray, params: Tuple[float]) -> np.ndarray:
        agentsNum = positionX.shape[0]
        boundaryLength, halfBoundaryLength, distanceD0 = params

        Pi = np.zeros(agentsNum)
        for i in range(agentsNum):
            xDiff = np.abs(positionX[:, 0] - positionX[i, 0])
            yDiff = np.abs(positionX[:, 1] - positionX[i, 1])
            neighborIdxs = np.where(
                (xDiff < distanceD0) | (boundaryLength - xDiff < distanceD0) & 
                (yDiff < distanceD0) | (boundaryLength - yDiff < distanceD0)
            )[0]
            if neighborIdxs.size == 0:
                continue

            subX = positionX[i] - positionX[neighborIdxs]
            deltaX = positionX[i] - (
                positionX[neighborIdxs] * (-halfBoundaryLength <= subX) * (subX <= halfBoundaryLength) + 
                (positionX[neighborIdxs] - boundaryLength) * (subX < -halfBoundaryLength) + 
                (positionX[neighborIdxs] + boundaryLength) * (subX > halfBoundaryLength)
            )
            distance = np.sqrt(np.sum(deltaX**2, axis=1))
            A = np.where((distance <= distanceD0) & (distance > 0))[0]

            Pi[i] = np.sum(1 / distance[A]) / (2 * np.pi)

        return Pi

    @property
    def rhoI(self):
        return self._neighbor_count_nI(self.positionX, self.dotRParams) / self.sensingArea

    @property
    def partiDiffD(self):
        if self.sensingStyle == "neighbor_count":
            return np.where(self.rhoI > self.rhoThres, self.diffusionDmax, self.diffusionDmin)
        elif self.sensingStyle == "inverse_distance_sum":
            return np.where(self._Pi(self.positionX, self.dotRParams) > self.Pth, self.diffusionDmax, self.diffusionDmin)

    @property
    def dotR(self):
        # Return a velocity-like white-noise term. update() multiplies dotR by dt,
        # so this gives the Euler-Maruyama displacement sqrt(2 D dt) * eta.
        return (
            np.sqrt(2 * self.partiDiffD[:, None] / self.dt) *
            np.random.standard_normal(self.positionX.shape)
        )

    def update(self):
        self.temp["dotR"] = self.dotR

        self.positionX = np.mod(self.positionX + self.temp["dotR"] * self.dt, self.boundaryLength)

    def append(self):
        if self.store is not None:
            if self.counts % self.shotsnaps != 0:
                return
            self.store.append(key="positionX", value=pd.DataFrame(self.positionX))

    def plot(self, ax: plt.Axes = None, particleStyle: str = 'scatter'):

        assert particleStyle in ['quiver', 'scatter']

        if ax is None:
            _, ax = plt.subplots(figsize=(5, 5))

        if particleStyle == 'scatter':
            ax.scatter(self.positionX[:, 0], self.positionX[:, 1], s=3, color="black", linewidths=0)
        elif particleStyle == 'quiver':
            norm = np.linalg.norm(self.temp["dotR"], axis=1, keepdims=True)
            norm[norm == 0] = 1
            unitDirection = self.temp["dotR"] / norm
            
            ax.quiver(
                self.positionX[:, 0], self.positionX[:, 1],
                unitDirection[:, 0], unitDirection[:, 1],
                angles="xy", scale_units="xy"
            )
        ax.set_xlim(0, self.boundaryLength)
        ax.set_ylim(0, self.boundaryLength)

    def state_description(self):
        return (
            pd.Series(self.partiDiffD).value_counts()
            .sort_index().rename(f"thresLp={self.thresLp:.2f}").to_frame()
        )

    def __str__(self):
        return (
            f"QuorumSensing_distanceDc{self.distanceDc:.2f},"
            f"thresLp{self.thresLp:.2f},boundaryLength{self.boundaryLength},"
            f"diffusionDmax{self.diffusionDmax},diffusionDmin{self.diffusionDmin},"
            f"style_{self.sensingStyle},"
            f"agentsNum{self.agentsNum},dt{self.dt},s{self.shotsnaps}"
            f"seed{self.randomSeed}"
        )