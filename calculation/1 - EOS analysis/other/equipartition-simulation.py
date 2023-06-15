# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes.constant import CALCULATION_DIR
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from copy import deepcopy
from tqdm import tqdm
import numpy as np
import os


# %%-------------------------------------   FUNCTION DEFINITION                 -------------------------------------> #
class system:

    def __init__(self, n_points=10, overall_k=15.):

        self.t = 0
        self.n_points = n_points
        self.point_pos = np.linspace(0, 1, num=(self.n_points + 2))

        self.point_speed = np.zeros(self.n_points + 2)

        self.__init_interaction_matrices()

        max_i = np.argmax(self.point_mass[1:-1])
        self.point_speed[max_i + 1] = overall_k / self.point_mass[max_i + 1]

        self.timeline = {

            "t": [self.t],
            "pos": [deepcopy(self.point_pos[1:-1])],
            "spd": [deepcopy(self.point_speed[1:-1])]

        }

    def __init_interaction_matrices(self):

        self.int_matrices = list()

        self.point_mass = np.random.rand(self.n_points + 2)
        self.point_mass[0] = np.inf
        self.point_mass[-1] = np.inf

        for i in range(1, self.n_points):

            m_1 = self.point_mass[i]
            m_2 = self.point_mass[i + 1]
            m_tot = m_1 + m_2
            dm = m_1 - m_2

            self.int_matrices.append(

                np.array([

                    [dm, 2*m_2],
                    [2*m_1, -dm]

                ]) / m_tot

            )

    def step(self):

        distances = self.point_pos[1:] - self.point_pos[:-1]
        rel_speed = self.point_speed[:-1] - self.point_speed[1:]
        inv_t_contacts = rel_speed / distances

        max_i = np.argmax(inv_t_contacts)
        dt = 1 / inv_t_contacts[max_i]

        self.t += dt
        self.point_pos += self.point_speed * dt

        if max_i == 0:

            self.point_pos[max_i + 1] = 0
            self.point_speed[max_i + 1] = -self.point_speed[max_i + 1]

        elif max_i == self.n_points:

            self.point_pos[max_i] = 1
            self.point_speed[max_i] = -self.point_speed[max_i]

        else:

            self.point_pos[max_i] = self.point_pos[max_i + 1]
            self.point_speed[max_i: max_i + 2] = np.dot(

                self.int_matrices[max_i - 1],
                self.point_speed[max_i: max_i + 2]

            )

    def generate_animation(self, steps=200, frames=1000, filename="equipartition_simulation.gif"):

        filepath = os.path.join(CALCULATION_DIR, "1 - EOS analysis", "other", "output", filename)
        pbar = tqdm(desc="Calculating Steps", total=steps)
        self.timeline = {

            "t": [self.t],
            "pos": [deepcopy(self.point_pos[1:-1])],
            "spd": [deepcopy(self.point_speed[1:-1])]

        }

        for i in range(steps):

            self.step()
            self.timeline["t"].append(self.t)
            self.timeline["pos"].append(deepcopy(self.point_pos[1:-1]))
            self.timeline["spd"].append(deepcopy(self.point_speed[1:-1]))

            pbar.update(1)

        pbar.close()

        fig, ax, sc = self.__init_plot()
        dt_tot = self.timeline["t"][-1] - self.timeline["t"][0]
        pbar = tqdm(desc="Generating Frames", total=(frames + 3))

        def update(frame):

            time = frame / (frames + 1) * dt_tot + self.timeline["t"][0]
            pbar.update(1)
            return self.__update_scatter(time),

        intr = 1 / (frames / 100)
        anim = animation.FuncAnimation(fig, update, frames=frames, interval=intr, blit=True)
        anim.save(filepath)
        pbar.close()

    def __init_plot(self):

        fig, ax = plt.subplots()
        self.sc = ax.scatter(

            self.point_pos[1:-1],
            np.zeros(self.n_points),
            s=self.point_mass[1:-1] * 80 + 10,
            alpha=0.5)

        ax.set_xlim([0, 1])
        return fig, ax, self.sc

    def __update_scatter(self, time):

        positions = self.point_pos[1:-1]

        if len(self.timeline["t"]) > 1:

            if self.timeline["t"][0] < time < self.timeline["t"][-1]:

                for i in range(len(self.timeline["t"])):

                    if time < self.timeline["t"][i]:

                        t_0 = self.timeline["t"][i - 1]
                        t_1 = self.timeline["t"][i]
                        x = (time - t_0) / (t_1 - t_0)

                        p_0 = self.timeline["pos"][i - 1]
                        p_1 = self.timeline["pos"][i]

                        positions = p_1 * x + p_0 * (1 - x)

                        break

        new_points = np.array([

            positions,
            np.zeros(self.n_points)

        ]).T

        self.sc.set_offsets(new_points)
        return self.sc


# %%-------------------------------------   GENERATE ANIMATION                  -------------------------------------> #
n_points = 25
k = 30 * (n_points / 10)
sys = system(n_points=n_points, overall_k=k)
sys.generate_animation(steps=1500, frames=2000)
