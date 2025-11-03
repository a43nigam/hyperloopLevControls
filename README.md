# hyperloopLevControls

MATLAB code for levitation controls for Cornell Hyperloop's minipod. The one_mag controls are for a 1 DoF (z) system with one magnet. The 4 mag controls are for a 3 DoF system with 4 magnets. 

Both codes use PID, the coefficients are hardcoded and work decently well, I didn't really fine tune them. 

The 4 mag code assumes that magnets are placed at 4 corners of some chassis, with chassis parameters editable. The idea is that we only need to know the pitch, yaw, and z of the center of mass of the chassis; assuming the chassis and all fixtures (like magnets) are rigid, we can easily use those parameters to figure out z displacement of each of the magnets. This doesn't take into account the tilt of magnets, but we assume that by the time tilt is not negligable, the pitch or yaw will be so great that at least one of the magnets will have struck the track, making it pointless to try to account for any non-vertical (shear) forces.

Note that both codes include an inductance factor - coils don't charge or discharge immediately, they act as inductors. We account for this by estimating their resistance (R) and inductance (L). The 4 magnet controls also have a time step factor, which discretizes current updates provided to the motor. The ODE solver automatically updates all values nearly continuously, but this is not accurate to real life, where we update our current in discrete intervals.

In both codes, we add disturbances at different time intervals to see how the system reacts

<img width="1556" height="1002" alt="Screenshot 2025-11-03 at 10 41 14 AM" src="https://github.com/user-attachments/assets/a4102f21-ef14-4cca-a4a1-ba0460ed3534" />
One magnet controls

<img width="1520" height="998" alt="Screenshot 2025-11-03 at 10 41 28 AM" src="https://github.com/user-attachments/assets/954b30c1-0e72-4281-8d40-1a124cc8cc60" />
4 magnet controls
