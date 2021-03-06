import tkinter
import time
from PIL import Image
import json
import ctypes

# width of the animation window
import numpy as np
from physics_system import PhysicsSystem
import errno
import os
from datetime import datetime

animation_window_width = 800
# height of the animation window
animation_window_height = 600
# initial x position of the ball
animation_ball_start_xpos = 50
# initial y position of the ball
animation_ball_start_ypos = 50
# radius of the ball
animation_ball_radius = 30
# the pixel movement of ball for each iteration
animation_ball_min_movement = 5
# delay between successive frames in seconds
animation_refresh_seconds = 0.01


# The main window of the animation
def create_simulation_window(bg_color):

    window = tkinter.Tk()
    window.title("Physics Simulation")
    # Uses python 3.6+ string interpolation
    window.geometry(f'{animation_window_width}x{animation_window_height}')
    window.configure(bg=bg_color)
    return window


# Create a canvas for animation and add it to main window
def create_animation_canvas(window, bg_color):
    canvas = tkinter.Canvas(window)
    canvas.configure(bg=bg_color)
    canvas.pack(fill="both", expand=True)  # this tells the canvas to take up the entire window
    return canvas


class PhysicsEngine:
    def __init__(self, physics_system):
        self.physics_system = physics_system
        self.background_color = "black"
        self.window_size = [600, 600]
        self.window = create_simulation_window("pink")
        self.canvas_size = [500, 500]
        self.canvas = create_animation_canvas(self.window, "white")
        self.origin = [0, 0]
        self.scale = 1

    def run_interval(self, delta_t=0.1, stop_time=10, sample_times=None):
        t = self.physics_system.global_time

        sample_times.sort()
        samples = []
        while t < stop_time:
            while len(sample_times) > 0 and sample_times[0] < t + (delta_t / 2):
                sample_time = sample_times.pop(0)
                if abs(t - sample_time) <= (delta_t / 2):
                    particles = [particle.get_object() for particle in self.physics_system.particles]
                    sample = {"time": sample_time, "particles": particles}
                    samples.append(sample)

            self.physics_system.time_step(delta_t)
            t = self.physics_system.global_time

        if sample_times is not None:
            timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
            f = open(f"./samples/{timestamp}.json", "w")
            json.dump(samples, f, indent=2)
            f.close()
            return timestamp

    def live_animation(self, framerate=60, updates_per_frame=10, sample_times=None):
        # The actual execution starts here
        particles = self.physics_system.particles
        particle_to_sprite = {}

        # creates a sprite for each of  the particles in the system
        for particle in particles:
            x = particle.position[0]
            y = particle.position[1]
            r = particle.radius
            y = 600 - y  # this corrects for the fact that the canvas y-axis is inverted and originates from the top
            particle_to_sprite[particle] = self.canvas.create_oval(x - r, y - r, x + r, y + r, fill="black")

        while True:
            prev_positions = {}
            for particle in particles:
                prev_positions[particle] = particle.position.copy()

            for _ in range(updates_per_frame):
                self.physics_system.time_step(1 / (framerate * updates_per_frame))

            position_deltas = {}
            for particle in particles:
                position_deltas[particle] = particle.position - prev_positions[particle]

            for particle in particles:
                sprite = particle_to_sprite[particle]
                position_delta = position_deltas[particle]
                self.canvas.move(sprite, position_delta[0], -1 * position_delta[1])

            self.window.update()

            if sample_times is not None:
                while len(sample_times) > 0 and (sample_times[0] <= self.physics_system.global_time):
                    t = sample_times.pop(0)
                    if t >= self.physics_system.global_time - (1 / framerate):
                        self.canvas.update()
                        filename = "../samples/" + str(t)
                        canvas_file = filename + ".ps"
                        self.canvas.postscript(file=canvas_file, colormode='color')
                        img = Image.open(canvas_file)
                        image_file = filename + ".png"
                        img.save(image_file)

            time.sleep(1 / framerate)
