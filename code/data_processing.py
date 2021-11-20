import json
import tkinter
from PIL import Image
import os

def json_to_images(timestamp):

    f = open(f"samples/{timestamp}.json", "r")
    samples = json.load(f)

    if not os.path.exists(f"samples/{timestamp}"):
        os.mkdir(f"samples/{timestamp}")

    for sample in samples:
        # get the data from the JSON object
        time = sample["time"]
        particles = sample["particles"]

        # create the canvas for drawing the system
        canvas = tkinter.Canvas(width=128, height=128)
        canvas.configure(bg="white")

        # draw the particles on the canvas
        for particle in particles:
            x = particle["position"][0]
            y = particle["position"][1]
            r = particle["radius"]
            canvas.create_oval(x - r, y - r, x + r, y + r, fill="black")

        canvas.update()
        filename = f"samples/{timestamp}/{str(time)}"
        canvas.postscript(file=f"{filename}.ps", colormode='color')
        img = Image.open(f"{filename}.ps")
        img.save(f"{filename}.png")
        os.remove(f"{filename}.ps")
