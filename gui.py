# from tkinter import *
# root = Tk()
#
# # the_label = Label(root, text="This is too easy")
# # the_label.pack()
#
# top_frame = Frame(root)
# top_frame.pack()
#
# bottom_frame = Frame(root)
# bottom_frame.pack(side=BOTTOM)
#
# button1 = Button(top_frame, text="Button 1", fg="red")
# button2 = Button(top_frame, text="Button 2", fg="green")
# button3 = Button(top_frame, text="Button 3", fg="blue")
# button4 = Button(bottom_frame, text="Button 4", fg="black")
#
# button1.pack(side=LEFT)
# button2.pack(side=LEFT)
# button3.pack(side=LEFT)
# button4.pack()
#
# root.mainloop()



import numpy as np
import sys
if sys.version_info[0] < 3:
    import Tkinter as tk
else:
    import tkinter as tk
import matplotlib.backends.tkagg as tkagg
from matplotlib.backends.backend_agg import FigureCanvasAgg


def draw_figure(canvas, figure, loc=(0, 0)):
    """ Draw a matplotlib figure onto a Tk canvas

    loc: location of top-left corner of figure on canvas in pixels.
    Inspired by matplotlib source: lib/matplotlib/backends/backend_tkagg.py
    """
    figure_canvas_agg = FigureCanvasAgg(figure)
    figure_canvas_agg.draw()
    figure_x, figure_y, figure_w, figure_h = figure.bbox.bounds
    figure_w, figure_h = int(figure_w), int(figure_h)
    photo = tk.PhotoImage(master=canvas, width=figure_w, height=figure_h)

    # Position: convert from top-left anchor to center anchor
    canvas.create_image(loc[0] + figure_w/2, loc[1] + figure_h/2, image=photo)

    # Unfortunately, there's no accessor for the pointer to the native renderer
    tkagg.blit(photo, figure_canvas_agg.get_renderer()._renderer, colormode=2)

    # Return a handle which contains a reference to the photo object
    # which must be kept live or else the picture disappears
    return photo


def init_draw(fig1, fig2):
    # Create a canvas
    w, h = 1600, 1600
    window = tk.Tk()
    window.title("Genetic Algorithm (GA) based structural search algorithms")
    canvas1 = tk.Canvas(window, width=w/2, height=h/2)
    canvas2 = tk.Canvas(window, width=w/2, height=h/2)
    canvas1.pack(side=tk.LEFT)
    canvas2.pack(side=tk.LEFT)


    # Add more elements to the canvas, potentially on top of the figure
    # canvas.create_line(200, 50, fig_x + fig_w / 2, fig_y + fig_h / 2)
    # canvas.create_text(200, 50, text="Zero-crossing", anchor="s")

    # Let Tk take over


    fig_photo_1 = draw_figure(canvas1, fig1)
    fig_photo_2 = draw_figure(canvas2, fig2)
    fig_w_1, fig_h_1 = fig_photo_1.width(), fig_photo_1.height()
    fig_w_2, fig_h_2 = fig_photo_2.width(), fig_photo_2.height()

    tk.mainloop()