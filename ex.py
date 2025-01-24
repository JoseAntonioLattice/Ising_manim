
import numpy as np
import matplotlib.pyplot as plt
from IPython import display
import time

pause_time = 0.2  # seconds between frames
a = np.random.rand(3,3)

for t in range(0,10):
    plt.imshow(a)
    plt.title(t)
    display.display(plt.gcf())
    display.clear_output(wait=True)
    time.sleep(pause_time)

    a = np.random.rand(3,3)
