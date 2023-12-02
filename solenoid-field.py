import matplotlib.pyplot as plt
import numpy as np

# Solenoid parameters
SOLENOID_HEIGHT = 5
SOLENOID_COILS = 48
SOLENOID_RADIUS = 1
SOLENOID_NUM_POINTS = 1500
SOLENOID_CURRENT = 0.02

# Plot parameters
PLOT_HORIZ_MARGIN = 0.5
PLOT_VERT_MARGIN = 3
PLOT_VERT_LIMIT = SOLENOID_HEIGHT + PLOT_VERT_MARGIN
PLOT_HORIZ_LIMIT = 2 * SOLENOID_RADIUS + PLOT_HORIZ_MARGIN

# Vector field parameters
FIELD_HORIZ_STEP = 0.5
FIELD_VERT_STEP = 0.5

# Set up the figure and plot
fig = plt.figure()
fig.set_figwidth(10)
fig.set_figheight(8)
ax = fig.add_subplot(projection='3d')

# Set up the axes
ax.axes.set_xlim3d(left=-PLOT_HORIZ_LIMIT, right=PLOT_HORIZ_LIMIT) 
ax.axes.set_ylim3d(bottom=-PLOT_HORIZ_LIMIT, top=PLOT_HORIZ_LIMIT) 
ax.axes.set_zlim3d(bottom=-PLOT_VERT_LIMIT, top=PLOT_VERT_LIMIT) 

# Parametric definition of a solenoid
def solenoid():
  theta = np.linspace(-2 * SOLENOID_COILS * np.pi, 2 * SOLENOID_COILS * np.pi, SOLENOID_NUM_POINTS)
  z = np.linspace(- SOLENOID_HEIGHT / 2, SOLENOID_HEIGHT / 2, SOLENOID_NUM_POINTS)
  x = SOLENOID_RADIUS * np.cos(theta)
  y = SOLENOID_RADIUS * np.sin(theta)
  return x, y, z

# Calculates the magnetic field at a given point
def fieldAtPoint(_r):
  # If the point is too close to the coil, it will have an enormous magnitude, so don't draw it
  distance_from_xy_origin = _r[0]**2 + _r[1]**2
  if distance_from_xy_origin < (SOLENOID_RADIUS * 1.2)**2 and distance_from_xy_origin > (SOLENOID_RADIUS * 0.8)**2:
    return np.array([0, 0, 0])

  r = np.array(_r) # Field point

  # Line integral along the solenoid, using the point divisions as dl
  [curve_x, curve_y, curve_z] = solenoid()
  _b = np.array([0, 0, 0])
  for i in range(SOLENOID_NUM_POINTS - 1):
    # Get the approximate dl as a segment of the curve. This won't be perfect,
    # but it's close enough for the purposes of numerical integration here,
    # assuming SOLENOID_NUM_POINTS is large enough.
    dlx = curve_x[i + 1] - curve_x[i]
    dly = curve_y[i + 1] - curve_y[i]
    dlz = curve_z[i + 1] - curve_z[i]
    dl = np.array([dlx, dly, dlz])

    # Represents the displacement from the source point to the field point
    sr = r - np.array([curve_x[i], curve_y[i], curve_z[i]])
    srmag = np.sqrt(sr.dot(sr)) # get the magnitude of the sr vector

    # Biot-Savart law, with constants removed
    db = ((SOLENOID_CURRENT * np.cross(dl, sr)) / (srmag**3))

    # Add to the running total
    _b = _b + db
  return _b

# Draws the solenoid wire shape
def drawSolenoid():
  x, y, z = solenoid()
  ax.plot(x, y, z, label='Solenoid', c='black')
  ax.legend()

# Draws the magnetic field vectors produced by the vector field
def drawField():
  # Set up a 3D grid of equally spaced points in each dimension
  x_points = np.arange(-PLOT_HORIZ_LIMIT, PLOT_HORIZ_LIMIT, FIELD_HORIZ_STEP)
  y_points = np.arange(-PLOT_HORIZ_LIMIT, PLOT_HORIZ_LIMIT, FIELD_HORIZ_STEP)
  z_points = np.arange(-PLOT_VERT_LIMIT, PLOT_VERT_LIMIT, FIELD_VERT_STEP)

  # Count the number of points in each dimension, for later use
  num_x_points = len(x_points)
  num_y_points = len(y_points)
  num_z_points = len(z_points)

  # Generate the mesh grid from the points made above
  grid_x, grid_y, grid_z = np.meshgrid(x_points, y_points, z_points, indexing='ij')

  # Initialize vectors for each grid point as (0, 0, 0) to start
  vectors_x = [[[float(0) for _ in b] for b in a] for a in grid_x]
  vectors_y = [[[float(0) for _ in b] for b in a] for a in grid_y]
  vectors_z = [[[float(0) for _ in b] for b in a] for a in grid_z]

  # Iterate through every point in the grid, and calculate its field
  for i in range(num_x_points):
    for j in range(num_y_points):
      for k in range(num_z_points):
        # Get the point's coordinates
        point_x = grid_x[i, j, k]
        point_y = grid_y[i, j, k]
        point_z = grid_z[i, j, k]
        point = [point_x, point_y, point_z]

        # Find the field at that point
        b = fieldAtPoint(point)

        # Record the field in the respective component matrices
        vectors_x[i][j][k] = b[0]
        vectors_y[i][j][k] = b[1]
        vectors_z[i][j][k] = b[2]

  # Draw the vector field as a collection of arrows
  ax.quiver(grid_x, grid_y, grid_z, vectors_x, vectors_y, vectors_z, arrow_length_ratio = 0.02)

# Draw the components
drawSolenoid()
drawField()

# Render the plot
plt.show()