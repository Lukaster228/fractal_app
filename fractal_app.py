#Программа для генерации фракталов. Салия Лука Мерабович КММО-01-23

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox

#Генерация Мандельброта
def mandelbrot(c, max_iter):
    z = c
    for n in range(max_iter):
        if abs(z) > 2:
            return n
        z = z*z + c
    return max_iter

def mandelbrot_set(xmin, xmax, ymin, ymax, width, height, max_iter):
    X = np.linspace(xmin, xmax, width)
    Y = np.linspace(ymin, ymax, height)
    C = X[:, None] + 1j*Y[None, :]
    M = np.zeros(C.shape, dtype=int)
    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            M[i, j] = mandelbrot(C[i, j], max_iter)
    return M

#Генерация Жюлиа
def julia(c, z, max_iter):
    for n in range(max_iter):
        if abs(z) > 2:
            return n
        z = z*z + c
    return max_iter

def julia_set(c, xmin, xmax, ymin, ymax, width, height, max_iter):
    X = np.linspace(xmin, xmax, width)
    Y = np.linspace(ymin, ymax, height)
    Z = X[:, None] + 1j*Y[None, :]
    J = np.zeros(Z.shape, dtype=int)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            J[i, j] = julia(c, Z[i, j], max_iter)
    return J

#Генерация H-Tree
def draw_h_tree(ax, order, x, y, size):
    if order == 0:
        return
    x0, x1 = x - size / 2, x + size / 2
    y0, y1 = y - size / 2, y + size / 2
    ax.plot([x0, x1], [y, y], 'b')
    ax.plot([x0, x0], [y0, y1], 'b')
    ax.plot([x1, x1], [y0, y1], 'b')
    draw_h_tree(ax, order - 1, x0, y0, size / 2)
    draw_h_tree(ax, order - 1, x0, y1, size / 2)
    draw_h_tree(ax, order - 1, x1, y0, size / 2)
    draw_h_tree(ax, order - 1, x1, y1, size / 2)

#Генерация Y-Tree
def draw_y_tree(ax, order, x, y, size, angle):
    if order == 0:
        return
    x1, y1 = x + size * np.cos(angle), y + size * np.sin(angle)
    ax.plot([x, x1], [y, y1], 'b')
    draw_y_tree(ax, order - 1, x1, y1, size * 0.7, angle + np.pi / 6)
    draw_y_tree(ax, order - 1, x1, y1, size * 0.7, angle - np.pi / 6)

#Генерации Кривой Дракона
def draw_dragon_curve(ax, order, x0, y0, x1, y1):
    def dragon_recursive(ax, order, p, q, left):
        if order == 0:
            ax.plot([p[0], q[0]], [p[1], q[1]], 'b')
        else:
            mid = ((p[0] + q[0]) / 2 + (q[1] - p[1]) / 2 * (1 if left else -1),
                   (p[1] + q[1]) / 2 - (q[0] - p[0]) / 2 * (1 if left else -1))
            dragon_recursive(ax, order - 1, p, mid, True)
            dragon_recursive(ax, order - 1, mid, q, False)
    dragon_recursive(ax, order, (x0, y0), (x1, y1), True)

#Генерация Кривой Коха
def draw_koch_curve(ax, order, p1, p2):
    if order == 0:
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'b')
    else:
        third = ((2*p1[0] + p2[0]) / 3, (2*p1[1] + p2[1]) / 3)
        two_third = ((p1[0] + 2*p2[0]) / 3, (p1[1] + 2*p2[1]) / 3)
        peak = ((third[0] + two_third[0]) / 2 - (two_third[1] - third[1]) * np.sqrt(3) / 2,
                (third[1] + two_third[1]) / 2 + (two_third[0] - third[0]) * np.sqrt(3) / 2)
        draw_koch_curve(ax, order - 1, p1, third)
        draw_koch_curve(ax, order - 1, third, peak)
        draw_koch_curve(ax, order - 1, peak, two_third)
        draw_koch_curve(ax, order - 1, two_third, p2)

# GUI для выбора и отображения фракталов
class FractalApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Fractal Generator")
        self.create_widgets()
        self.scale = 1.0
        self.offset_x = 0
        self.offset_y = 0
        self.drag_start_x = None
        self.drag_start_y = None
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        
    def create_widgets(self):
        self.fractal_type = tk.StringVar(value="Mandelbrot")
        self.max_iter = tk.IntVar(value=100)
        self.real_part = tk.DoubleVar(value=-0.7)
        self.imag_part = tk.DoubleVar(value=0.27015)
        
        controls_frame = ttk.Frame(self.root)
        controls_frame.pack(side=tk.LEFT, fill=tk.Y)
        
        fractal_label = ttk.Label(controls_frame, text="Fractal Type")
        fractal_label.pack(pady=5)
        
        fractal_combo = ttk.Combobox(controls_frame, textvariable=self.fractal_type, values=["Mandelbrot", "Julia", "H-Tree", "Y-Tree", "Dragon Curve", "Koch Curve"])
        fractal_combo.pack(pady=5)
        
        iter_label = ttk.Label(controls_frame, text="Max Iterations")
        iter_label.pack(pady=5)
        
        iter_spinbox = ttk.Spinbox(controls_frame, from_=10, to=1000, textvariable=self.max_iter)
        iter_spinbox.pack(pady=5)
        
        real_label = ttk.Label(controls_frame, text="Real Part (for Julia)")
        real_label.pack(pady=5)
        
        real_entry = ttk.Entry(controls_frame, textvariable=self.real_part)
        real_entry.pack(pady=5)
        
        imag_label = ttk.Label(controls_frame, text="Imaginary Part (for Julia)")
        imag_label.pack(pady=5)
        
        imag_entry = ttk.Entry(controls_frame, textvariable=self.imag_part)
        imag_entry.pack(pady=5)
        
        generate_button = ttk.Button(controls_frame, text="Generate", command=self.generate_fractal)
        generate_button.pack(pady=20)
        
        self.canvas_frame = ttk.Frame(self.root)
        self.canvas_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.canvas = None
        
    def generate_fractal(self):
        if self.canvas:
            self.canvas.get_tk_widget().destroy()

        max_iter_str = self.max_iter.get()

        if max_iter_str < 1:
            messagebox.showerror("Error", "Max iterations must be greater than or equal to 1")
            return

        fig, self.ax = plt.subplots()
        self.draw_current_fractal()
        
        self.canvas = FigureCanvasTkAgg(fig, master=self.canvas_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.canvas.mpl_connect("scroll_event", self.on_scroll)
        self.canvas.mpl_connect("button_press_event", self.on_button_press)
        self.canvas.mpl_connect("button_release_event", self.on_button_release)
        self.canvas.mpl_connect("motion_notify_event", self.on_mouse_move)
        
    def draw_current_fractal(self):
        self.ax.clear()
        if self.fractal_type.get() == "Mandelbrot":
            M = mandelbrot_set(-2.0, 1.0, -1.5, 1.5, 800, 800, self.max_iter.get())
            self.ax.imshow(M, cmap='hot', extent=[-2, 1, -1.5, 1.5])
        elif self.fractal_type.get() == "Julia":
            J = julia_set(complex(self.real_part.get(), self.imag_part.get()), -1.5, 1.5, -1.5, 1.5, 800, 800, self.max_iter.get())
            self.ax.imshow(J, cmap='hot', extent=[-1.5, 1.5, -1.5, 1.5])
        elif self.fractal_type.get() == "H-Tree":
            self.ax.set_aspect('equal')
            draw_h_tree(self.ax, self.max_iter.get(), 0, 0, 200)
        elif self.fractal_type.get() == "Y-Tree":
            self.ax.set_aspect('equal')
            draw_y_tree(self.ax, self.max_iter.get(), 0, -200, 200, np.pi / 2)
        elif self.fractal_type.get() == "Dragon Curve":
            self.ax.set_aspect('equal')
            draw_dragon_curve(self.ax, self.max_iter.get(), -100, 0, 100, 0)
        elif self.fractal_type.get() == "Koch Curve":
            self.ax.set_aspect('equal')
            draw_koch_curve(self.ax, self.max_iter.get(), (-200, 0), (200, 0))
        self.ax.set_xlim([-200, 200])
        self.ax.set_ylim([-200, 200])
        self.apply_transformations()
        
    def apply_transformations(self):
        self.ax.set_xlim(self.offset_x - 200/self.scale, self.offset_x + 200/self.scale)
        self.ax.set_ylim(self.offset_y - 200/self.scale, self.offset_y + 200/self.scale)
        
    def on_scroll(self, event):
        if event.button == 'up':
            self.scale *= 1.1
        elif event.button == 'down':
            self.scale /= 1.1
        self.apply_transformations()
        self.canvas.draw()
        
    def on_button_press(self, event):
        if event.button == 1:
            self.drag_start_x = event.x
            self.drag_start_y = event.y
            
    def on_button_release(self, event):
        if event.button == 1:
            self.drag_start_x = None
            self.drag_start_y = None
            
    def on_mouse_move(self, event):
        if self.drag_start_x is not None and self.drag_start_y is not None:
            dx = event.x - self.drag_start_x
            dy = event.y - self.drag_start_y
            self.offset_x -= dx / self.scale
            self.offset_y -= dy / self.scale
            self.drag_start_x = event.x
            self.drag_start_y = event.y
            self.apply_transformations()
            self.canvas.draw()

    def on_closing(self):
        if self.canvas:
            self.canvas.get_tk_widget().destroy()
            plt.close('all')  # Закрыть все окна графиков matplotlib
        self.root.destroy()

if __name__ == "__main__":  
    root = tk.Tk()
    app = FractalApp(root)
    root.mainloop()
    
