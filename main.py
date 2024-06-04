from vpython import *

# Scene
global temperature_readout, pressure_readout, KE_readout, PE_readout, total_energy_readout
WIDTH = 70
WALL_STIFFNESS = 50
running = False

# Atom Attributes: quantity, x/y position, velocity, acceleration
atoms = []
N = 300
x = [0] * N
y = [0] * N
vx = [0] * N
vy = [0] * N
ax = [0] * N
ay = [0] * N

# Calculations
total_PE = 0
total_KE = 0
current_KE = 0
avg_KE = 0
total_force_on_walls = 0
avg_force_on_walls = 0
temperature = 0
pressure = 0
total_energy = 0


def run_simulation():
    global avg_force_on_walls, current_KE, total_KE, avg_KE

    initialize_scene()
    initialize_atoms()

    while True:
        rate(1500)
        if running:
            for i in range(10):
                time_step()
                avg_force_on_walls += 1
            for i in range(N):  #
                atoms[i].pos = vec(x[i], y[i], 0)
            current_KE = calculate_kinetic_energy()
            total_KE += current_KE
            avg_KE += N
            calculate_data()
            update_data_text()


def initialize_scene():
    scene.background = vec(0, 0, 0)
    scene.width = 600
    scene.height = scene.width
    scene.center = vec(WIDTH / 2, WIDTH / 2, 0)
    scene.range = WIDTH / 2
    scene.fov = 0.01
    scene.autoscale = False
    scene.userspin = False
    scene.userzoom = False

    create_buttons()
    create_data_text()


def create_buttons():
    scene.append_to_caption("\n")
    scene.append_to_caption("\t")
    button(text="Start/Stop", bind=start_stop_simulation)
    scene.append_to_caption("\t")
    button(text="Add Energy to System", bind=add_energy)
    scene.append_to_caption("\t")
    button(text="Remove Energy from System", bind=remove_energy)
    scene.append_to_caption("\n\n")


def create_data_text():
    global temperature_readout, pressure_readout, KE_readout, PE_readout, total_energy_readout

    temperature_readout = wtext(text="Average Temperature: 0")
    scene.append_to_caption("\n")

    pressure_readout = wtext(text="Average Pressure: 0")
    scene.append_to_caption("\n")

    KE_readout = wtext(text="System Kinetic Energy: {:.2f}".format(calculate_kinetic_energy()))
    scene.append_to_caption("\n")

    PE_readout = wtext(text="System Potential Energy: {:.2f}".format(total_PE))
    scene.append_to_caption("\n")

    total_energy_readout = wtext(text="Total System Energy: {:.2f}".format(0.0))


def initialize_atoms():
    y_position = WIDTH - 1
    x_position = 1

    for i in range(N):
        if x_position + 1 > WIDTH:
            x_position = 1
            y_position -= 1.1

        x[i] = x_position
        y[i] = y_position
        vx[i] = 0
        vy[i] = 0
        ax[i] = 0
        ay[i] = 0

        x_position += 1.1
        atoms.append(sphere(radius=0.5, pos=vec(x[i], y[i], 0), color=color.green))


def move_atoms():
    global total_PE, total_force_on_walls
    total_PE = 0

    for i in range(N):
        if x[i] < 0.5:
            bounce_left(i)
        elif x[i] > (WIDTH - 0.5):
            bounce_right(i)
        else:
            ax[i] = 0

        if y[i] < 0.5:
            bounce_up(i)
        elif y[i] > (WIDTH - 0.5):
            bounce_down(i)
        else:
            ay[i] = 0

        calculate_accelerations(i)


def calculate_accelerations(i):
    global total_PE, total_force_on_walls
    force_cutoff_squared = 9
    energy_offset = abs(4 * (1 / force_cutoff_squared ** 6 - 1 / force_cutoff_squared ** 3))

    for j in range(i):
        x_difference = x[i] - x[j]
        y_difference = y[i] - y[j]
        r2 = x_difference**2 + y_difference**2

        if r2 < force_cutoff_squared:
            force_factor = 24 * (2 * (1 / r2**7) - (1 / r2**4))
            x_force_factor = force_factor * x_difference
            y_force_factor = force_factor * y_difference

            ax[i] += x_force_factor
            ay[i] += y_force_factor
            ax[j] += -x_force_factor
            ay[j] += -y_force_factor

            total_PE += 4 * ((1 / r2**2) * (1 / r2**4) - (1 / r2**2) * (1/r2)) + energy_offset
        else:
            pass


def bounce_left(i):
    global total_PE, total_force_on_walls

    x_displacement = 0.5 - x[i]
    spring_force = WALL_STIFFNESS * x_displacement
    ax[i] = spring_force

    total_PE += 0.5 * spring_force * x_displacement
    total_force_on_walls += spring_force


def bounce_right(i):
    global total_PE, total_force_on_walls

    x_displacement = WIDTH - 0.5 - x[i]
    spring_force = WALL_STIFFNESS * x_displacement
    ax[i] = spring_force

    total_PE += 0.5 * spring_force * x_displacement
    total_force_on_walls += -spring_force


def bounce_up(i):
    global total_PE, total_force_on_walls

    y_displacement = 0.5 - y[i]
    spring_force = WALL_STIFFNESS * y_displacement
    ay[i] = spring_force

    total_PE += 0.5 * spring_force * y_displacement
    total_force_on_walls += spring_force


def bounce_down(i):
    global total_PE, total_force_on_walls

    y_displacement = WIDTH - 0.5 - y[i]
    spring_force = WALL_STIFFNESS * y_displacement
    ay[i] = spring_force

    total_PE += 0.5 * spring_force * y_displacement
    total_force_on_walls += -spring_force


def time_step():
    dt = 0.02

    for i in range(N):
        # Verlet Algorithm
        x[i] += vx[i] * dt + ax[i] * 0.5 * dt**2
        y[i] += vy[i] * dt + ay[i] * 0.5 * dt**2
        vx[i] += ax[i] * 0.5 * dt
        vy[i] += ay[i] * 0.5 * dt

    move_atoms()

    for i in range(N):
        vx[i] += ax[i] * 0.5 * dt
        vy[i] += ay[i] * 0.5 * dt


def start_stop_simulation():
    global running
    running = not running


def add_energy():
    reset_averages()
    for i in range(N):
        vx[i] *= 1.1
        vy[i] *= 1.1


def remove_energy():
    reset_averages()
    for i in range(N):
        vx[i] /= 1.1
        vy[i] /= 1.1


def reset_averages():
    global total_KE, avg_KE, total_force_on_walls, avg_force_on_walls
    total_KE = 0
    avg_KE = 0
    total_force_on_walls = 0
    avg_force_on_walls = 0


def calculate_kinetic_energy():
    global current_KE
    current_KE = 0

    for i in range(N):
        speed_squared = vx[i]**2 + vy[i]**2
        current_KE += 0.5 * speed_squared

    return current_KE


def calculate_data():
    global temperature, pressure, total_energy, current_KE
    temperature = total_KE / avg_KE

    if avg_force_on_walls == 0:
        pressure = 0
    else:
        pressure = total_force_on_walls / (avg_force_on_walls * WIDTH * 4)

    total_energy = total_PE + current_KE


def update_data_text():
    KE_readout.text = "System Kinetic Energy: {:.2f}".format(current_KE)
    PE_readout.text = "System Potential Energy: {:.2f}".format(total_PE)
    total_energy_readout.text = "Total System Energy: {:.2f}".format(total_energy)
    temperature_readout.text = "Average Temperature: {:.5f}".format(temperature)
    pressure_readout.text = "Average Pressure: {:.5f}".format(pressure)


run_simulation()
