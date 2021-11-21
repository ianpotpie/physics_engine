import json
import tkinter
from PIL import Image
import os


def data_to_images(timestamp, width=128, height=128):
    f = open(f"samples/{timestamp}.json", "r")
    samples = json.load(f)

    if not os.path.exists(f"samples/{timestamp}"):
        os.mkdir(f"samples/{timestamp}")

    window = tkinter.Tk()
    for sample in samples:
        # get the data from the JSON object
        time = sample["time"]
        particles = sample["particles"]

        # create the canvas for drawing the system
        canvas = tkinter.Canvas(window, width=width, height=height)
        canvas.pack(fill="both", expand=True)
        canvas.configure(bg="white")

        # draw the particles on the canvas
        for particle in particles:
            x = particle["position"][0]
            y = particle["position"][1]
            r = particle["radius"]
            y = height - y  # this corrects for the fact that the canvas y-axis is inverted and originates from the top
            canvas.create_oval(x - r, y - r, x + r, y + r, fill="black")

        window.update()
        canvas.update()

        # canvas.update()
        filename = f"samples/{timestamp}/{str(time)}"
        canvas.postscript(file=f"{filename}.ps", colormode='color')
        canvas.destroy()
        img = Image.open(f"{filename}.ps")
        img.save(f"{filename}.png")
        os.remove(f"{filename}.ps")

    window.destroy()
