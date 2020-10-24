import random
import math
from tkinter import *
import tkinter.font as font


def on_segment(p, q, r):
    """
    Given three colinear points p, q, r, check if point q lies on line segment 'pr'.
    """
    return min(p[0], r[0]) <= q[0] <= max(p[0], r[0]) and min(p[1], r[1]) <= q[1] <= max(p[1], r[1])


def ccw(a, b, c):
    """
    Returns 1 if points A,B,C are CCW, -1 if CW, and 0 if collinear
    """
    orient = (c[1] - a[1]) * (b[0] - a[0]) - (b[1] - a[1]) * (c[0] - a[0])

    if orient > 0:
        return 1
    elif orient < 0:
        return -1
    return 0


def distance(p1, p2):
    return ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2) ** 0.5


def cost(p1, p2):
    return distance(p1, p2)


def move_towards(p_i, p_f):  # O(1)
    line = Line(p_i, p_f)
    if line.norm() <= STEP_SIZE:
        return p_f
    proportion = STEP_SIZE / line.norm()
    return point_along_at(line, proportion)


def point_along_at(line, proportion):
    p_new = tuple()
    for i in [0,1]:
        p_new += (line.p1[i] + proportion * (line.p2[i] - line.p1[i]),)
    return p_new


class Line:
    def __init__(self, p1, p2):  # p1, p2 are each a tuple (p.x, p.y)
        self.__p1 = p1
        self.__p2 = p2

    @property
    def p1(self):
        return self.__p1

    @property
    def p2(self):
        return self.__p2

    def norm(self):
        return distance(self.p1, self.p2)

    def intersects(self, l2):  # O(1)
        """
        Return true if and only if line segments AB and CD intersect
        """
        a, b = self.p1, self.p2
        c, d = l2.p1, l2.p2
        o1, o2, o3, o4 = ccw(a, b, c), ccw(a, b, d), ccw(a, c, d), ccw(b, c, d)

        # General Case
        if o1 != o2 and o3 != o4:
            return True

        # Special Case
        if (not o1 and on_segment(a, c, b)) or (not o2 and on_segment(a, d, b)) \
                or (not o3 and on_segment(c, a, d)) or (not o4 and on_segment(c, b, d)):
            return True
        return False


class Polygon:
    def __init__(self, points):
        self.__points = points
        self.__lines = []
        points_offset = points[1:] + [points[0]]
        for pair_points in zip(points, points_offset):
            self.lines.append(Line(*pair_points))

    @property
    def points(self):
        return self.__points

    @property
    def lines(self):
        return self.__lines

    def intersects(self, line):   # O(num vertices in polygon)
        for l2 in self.lines:
            if line.intersects(l2):
                return True
        return False

    def contains(self, point, size):   # O(num vertices in polygon)
        """
        Returns True if point p1 is inside the polygon.
        *Time complexity: O(n) where n is the vertices in polygon
        """
        if len(self.lines) < 3:
            return False

        extended_line = Line(point, (size[0]+1, point[1]))
        count = 0
        for line in self.lines:
            if extended_line.intersects(line) and ccw(point, line.p1, line.p2):
                count += 1
        return count % 2 == 1


class Rectangle(Polygon):
    def __init__(self, x1, y1, x2, y2):
        points = [(x1, y1),(x1, y2),(x2, y2),(x2, y1)]
        Polygon.__init__(self, points)


class Boundaries:
    def __init__(self, size, polygons):
        self.__size = size
        self.__polygons = polygons

    @property
    def size(self):
        return self.__size

    @property
    def polygons(self):
        return self.__polygons

    def intersects(self, line):  # O(num vertices in boundaries)
        for polygon in self.polygons:
            if polygon.intersects(line):
                return True
        return False

    def contains(self, point):   # O(num vertices in boundaries)
        for polygon in self.polygons:
            if polygon.contains(point, self.size):
                return True
        return False


def get_nearest(vertices, goal_point):   # O(num vertices in V)
    min = TOTAL_SIZE*2
    min_v = (-1,-1)
    for v in vertices:
        d = distance(v, goal_point)
        if d < min:
            min = d
            min_v = v
    return min_v


def get_near(vertices, goal_point, dist):   # O(num vertices in V)
    return [v for v in vertices if distance(v, goal_point) < dist]


PI = 3.14
TOTAL_SIZE = 1000
STEP_SIZE = 30  # eta
NUM_DIM = 2  # d
MU_X_FREE = TOTAL_SIZE ** NUM_DIM # Lebesgue measure of X_free space, upper bounded by X_sample
ZETA_DIM = PI * 1 * 1 #area of ball in d dimensions
GAMMA_RRT_STAR = 2 * ((1 + 1/NUM_DIM) * MU_X_FREE / ZETA_DIM)**(1/NUM_DIM)

# Visualization for both RRT and RRT*
class Visualize_Path_Finding(Frame):

    def __init__(self, boundaries, start, goal):
        super().__init__()

        self.start = start
        self.goal = goal
        self.boundaries = boundaries
        self.algorithm = 'RRT' #'RRT' or 'RRT*'
        self.init_vars()

        self.cost_frame = Frame(self)
        self.cost_frame.pack(side="top")
        self.canvas = Canvas(self)
        self.canvas.pack(fill=BOTH, expand=1)

        self.initUI()

    def init_vars(self):
        self.V_adj = {start: set()}
        self.parents = {start: start}
        self.costs = {start: 0}
        self.p_best_in_target = None

        self.visual_V_adj = {start: {}}
        self.prev_target = [None] * 2
        self.prev_p_random = None
        self.prev_edge = None
        self.cur_total_cost = None

    def initUI(self):

        self.master.title("Path Planning Algorithms (RRT and RRT*)")
        self.pack(fill=BOTH, expand=1)

        self.initBoundaries()
        self.createWidgets()
        self.show_vertex(start, 'red')
        self.show_goal(15)

    def initBoundaries(self):
        for polygon in self.boundaries.polygons:
            points_flattened = []
            for point in polygon.points:
                points_flattened.append(point[0])
                points_flattened.append(point[1])
            self.canvas.create_polygon(points_flattened, outline='blue', fill='blue', width=2)

    def createWidgets(self):
        myFont = font.Font(size=30)

        # Algorithm buttons
        self.RRT_frame = Frame(self)
        self.RRT_frame.pack(side="left")
        button_text = Label(self.RRT_frame, text="    Algorithm:", font=(None,30))
        button_text.pack(side=LEFT)
        button = [None] * 2
        def make_alg_button(i):
            text = 'RRT*' if i else 'RRT'
            button[i] = Button(self.RRT_frame, text=text, command=lambda: self.change_to_algorithm(text))
            button[i].pack(side = LEFT)
            button[i]['font'] = myFont
        for i in range(2):
            make_alg_button(i)

        # Node buttons
        self.widget_frame = Frame(self)
        self.widget_frame.pack(side="left")
        button_text = Label(self.widget_frame, text="    Number of Nodes to Add:", font=(None,30))
        button_text.pack(side=LEFT)
        button = [None] * 4
        def make_node_button(i):
            button[i] = Button(self.widget_frame, text=("1" + "0"*i),
                               command= lambda: self.add_vertices_RRT_Star(10**i) if self.algorithm == 'RRT*' else self.add_vertices_RRT(10**i))
            button[i].pack(side = LEFT)
            button[i]['font'] = myFont
        for i in range(4):
            make_node_button(i)

        # Exploration bias slider
        self.exp_frame = Frame(self)
        self.exp_frame.pack(side="left")
        self.exp_text = Label(self.exp_frame, text="    Exploration Bias:", font=(None, 30))
        self.exp_text.pack(side=LEFT)
        self.exploration_bias_slider = Scale(self.exp_frame, from_=0, to=1, resolution = 0.01, orient = HORIZONTAL, font=(None,30))
        self.exploration_bias_slider.pack(side = LEFT)

        # Goal radius slider
        self.goal_frame = Frame(self)
        self.goal_frame.pack(side="left")
        self.goal_text = Label(self.goal_frame, text="    Goal Radius:", font=(None, 30))
        self.goal_text.pack(side=LEFT)
        self.goal_radius_slider = Scale(self.goal_frame, from_=0, to=50, resolution=0.1, orient=HORIZONTAL,
                                             font=(None, 30), command = self.show_goal)
        self.goal_radius_slider.set(15)
        self.goal_radius_slider.pack(side=LEFT)

        # Cost text
        self.cost_text = Label(self.cost_frame, text="", font=(None, 30))
        self.update_cost_text()
        self.cost_text.pack(side=LEFT)

    def update_cost_text(self):
        new_text = self.algorithm + ': '
        new_text += str(len(self.V_adj)) + ' Node' + ('s' if len(self.V_adj) != 1 else '') + ', '
        new_text += ('Total Distance to Goal: ' + '{:.0f}'.format(self.cur_total_cost) + " pixels") if self.cur_total_cost else "Goal not Found"
        self.cost_text.config(text=new_text)

    def change_to_algorithm(self, new_alg):
        self.canvas.delete("vertices")
        self.canvas.delete("edges")
        self.init_vars()
        self.algorithm = new_alg

        self.update_cost_text()
        self.show_vertex(start, 'red')
        self.show_goal()

    def show_vertex(self, loc, color='firebrick4', size = 5):
        x1, y1 = (loc[0] - size), (loc[1] - size)
        x2, y2 = (loc[0] + size), (loc[1] + size)
        self.canvas.create_oval(x1, y1, x2, y2, fill=color, tags = ("vertices",))

    def show_target(self, loc, color='firebrick4', size = 5):
        x1, y1 = (loc[0] - size), (loc[1] - size)
        x2, y2 = (loc[0] + size), (loc[1] + size)
        self.prev_target[1] = self.prev_target[0]
        self.prev_target[0] = self.canvas.create_oval(x1, y1, x2, y2, fill=color, tags = ("vertices",))

    def show_once_vertex(self, loc, color='black', size = 5):
        x1, y1 = (loc[0] - size), (loc[1] - size)
        x2, y2 = (loc[0] + size), (loc[1] + size)
        self.prev_p_random = self.canvas.create_oval(x1, y1, x2, y2, fill=color, tags = ("vertices",))

    def show_goal(self, target_radius=None):
        if target_radius is None: target_radius = self.goal_radius_slider.get()
        for i in range(2):
            if self.prev_target[i]:
                self.canvas.delete(self.prev_target[i])
        self.show_target(self.goal, color = 'yellow', size=float(target_radius))
        self.show_target(self.goal, color = 'green')

    def show_edge(self, v1, v2, color = 'red2', width=1, tags = ("edges",)):
        x1, y1 = v1[0], v1[1]
        x2, y2 = v2[0], v2[1]
        self.visual_V_adj.setdefault(v1, {})[v2] = self.canvas.create_line(x1, y1, x2, y2, fill=color, width=width, tags=tags)
        self.visual_V_adj.setdefault(v2, {})[v1] = self.visual_V_adj[v1][v2]

    def show_once_edge(self, v1, v2, color = 'black', width=1):
        x1, y1 = v1[0], v1[1]
        x2, y2 = v2[0], v2[1]
        self.prev_edge = self.canvas.create_line(x1, y1, x2, y2, dash = (4,2), fill=color, width=width, tags = ("edges",))

    def show_path(self, path):
        for edge in zip(path, path[1:]):
            self.show_edge(*edge, 'lime green', width = 5, tags = ("goal edges","edges"))

    def in_target(self, node):
        return self.near(node, self.goal_radius_slider.get())

    def near(self, node, dist):
        return distance(node, self.goal) <= dist

    def path_to_goal(self, goal_node):
        # Construct path in reverse using parents dict
        self.canvas.delete("goal edges")
        path = []
        cur = goal_node
        while cur is not self.start:
            path.append(cur)
            cur = self.parents[cur]
        path.append(self.start)
        path.reverse()
        self.show_path(path)

    def add_vertices_RRT(self, num_times=1):   # O(num nodes * num_times)
        for count in range(num_times):
            # Choose random point not in boundaries
            if self.prev_p_random is not None:
                self.canvas.delete(self.prev_p_random)
                self.canvas.delete(self.prev_edge)
            while True:
                p_random = (random.uniform(0, self.boundaries.size[0]),
                            random.uniform(0, self.boundaries.size[1]))
                if not self.boundaries.contains(p_random):
                    break
            if random.uniform(0,1) < self.exploration_bias_slider.get():
                p_random = self.goal
            self.show_once_vertex(p_random)

            # Create new point by moving towards the chosen random point
            p_nearest = get_nearest(self.V_adj, p_random)
            p_new = move_towards(p_nearest, p_random)
            if p_new not in self.V_adj and not self.boundaries.intersects(Line(p_nearest, p_new)):
                self.V_adj.setdefault(p_nearest, set()).add(p_new)
                self.V_adj.setdefault(p_new, set()).add(p_nearest)
                self.parents[p_new] = p_nearest
                self.costs[p_new] = self.costs[p_nearest] + cost(p_new, p_nearest)

                self.show_vertex(p_new)
                self.show_edge(p_nearest, p_new)
                self.show_once_edge(p_new, p_random)
                if self.in_target(p_new) and (self.p_best_in_target is None or self.costs[p_new] < self.costs[self.p_best_in_target]):
                    self.p_best_in_target = p_new

        # Update goal path and cost text
        self.canvas.delete("goal edges")
        if self.p_best_in_target and self.in_target(self.p_best_in_target):
            self.cur_total_cost = self.costs[self.p_best_in_target]
            self.path_to_goal(self.p_best_in_target)
        else:
            self.p_best_in_target = None
            self.cur_total_cost = None
        self.update_cost_text()

    def update_costs(self, v):
        self.costs[v] = self.costs[self.parents[v]] + cost(self.parents[v], v)
        for neighbor in self.V_adj[v]:
            if neighbor != self.parents[v]:
                self.update_costs(neighbor)

    def add_vertices_RRT_Star(self, num_times=1):   # O(num nodes * num_times)
        for count in range(num_times):
            # Choose random point not in boundaries
            if self.prev_p_random is not None:
                self.canvas.delete(self.prev_p_random)
                self.canvas.delete(self.prev_edge)
            while True:
                p_random = (random.uniform(0, self.boundaries.size[0]),
                            random.uniform(0, self.boundaries.size[1]))
                if not self.boundaries.contains(p_random):
                    break
            if random.uniform(0,1) < self.exploration_bias_slider.get():
                p_random = self.goal
            self.show_once_vertex(p_random)

            # Create new point by moving towards the chosen random point
            p_nearest = get_nearest(self.V_adj, p_random)
            p_new = move_towards(p_nearest, p_random)
            if p_new not in self.V_adj and not self.boundaries.intersects(Line(p_nearest, p_new)):
                RRT_Star = min(GAMMA_RRT_STAR * (math.log(len(self.V_adj)) / len(self.V_adj))**(1/NUM_DIM), 2.5*STEP_SIZE)
                P_near = get_near(self.V_adj, p_new, RRT_Star)

                # Find minimum cost to reach p_new
                p_min, c_min = p_nearest, self.costs[p_nearest] + cost(p_nearest, p_new)
                for p_near in P_near:
                    if not self.boundaries.intersects(Line(p_near, p_new))\
                            and self.costs[p_near] + cost(p_near, p_new) < c_min:
                        c_min = self.costs[p_near] + cost(p_near, p_new)
                        p_min = p_near
                self.V_adj.setdefault(p_min, set()).add(p_new)
                self.V_adj.setdefault(p_new, set()).add(p_min)
                self.parents[p_new] = p_min
                self.costs[p_new] = c_min

                self.show_vertex(p_new)
                self.show_edge(p_min, p_new)
                self.show_once_edge(p_new, p_random)

                # Rewire the tree with updated minimum costs through p_new
                for p_near in P_near:
                    if not self.boundaries.intersects(Line(p_new, p_near))\
                            and self.costs[p_new] + cost(p_new, p_near) < self.costs[p_near]:
                        self.V_adj[p_near].remove(self.parents[p_near])
                        self.V_adj[self.parents[p_near]].remove(p_near)
                        self.canvas.delete(self.visual_V_adj[p_near][self.parents[p_near]])
                        self.parents[p_near] = p_new
                        self.update_costs(p_near)
                        self.V_adj.setdefault(p_near, set()).add(p_new)
                        self.V_adj.setdefault(p_new, set()).add(p_near)
                        self.show_edge(p_near, p_new)

        # Update goal path and cost text
        self.p_best_in_target = None
        self.cur_total_cost = None
        self.canvas.delete("goal edges")
        for node in self.V_adj:
            if self.in_target(node) and \
                    (self.p_best_in_target is None or self.costs[node] <= self.costs[self.p_best_in_target]):
                self.p_best_in_target = node
        if self.p_best_in_target:
            self.cur_total_cost = self.costs[self.p_best_in_target]
            self.path_to_goal(self.p_best_in_target)
        self.update_cost_text()


# Create boundaries for visualization
shapes = []
shapes.append(Rectangle(350, 0, 400, 800))
shapes.append(Polygon([(600, 150), (550,200), (650,200)]))
shapes.append(Polygon([(1000,0), (800,350), (600,550), (850, 400)]))
shapes.append(Polygon([(500,500), (450,550), (500,700), (550, 550)]))
shapes.append(Rectangle(1100, 525, 1300, 550))
shapes.append(Rectangle(1100, 550, 1125, 350))
shapes.append(Rectangle(1100, 350, 1300, 375))
boundaries = Boundaries((1600,900),shapes)
start = (500,100)
goal = (1200,450)

root = Tk()
size = (1600, 900)
root.geometry('1600x900')
v_RRT = Visualize_Path_Finding(boundaries, start, goal)
root.mainloop()