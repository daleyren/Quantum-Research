from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
def tostr(x):
    return '{0:.1f}'.format(x)
import time

Nat = 2
a0 = basis(Nat,0)
a1 = basis(Nat,1)

def r_theta2(theta,phi):
    """
    Single qubit rotation around axis theta for angle phi
    result : qobj
    """
    return Qobj([[np.cos(phi / 2), -1j*np.sin(phi / 2)*(np.cos(theta)-1.0j*np.sin(theta))],
                     [-1.0j*np.sin(phi / 2)*(np.cos(theta)+1.0j*np.sin(theta)), np.cos(phi / 2)]])

def unitary_operator(op, densitymat):
    """
    apply operator to densitymat and return the new density matrix
    """
    return op*densitymat*op.dag() 

def ideal_measurement(densitymat):
  copy = np.array(densitymat)
  copy[0][1]=0
  copy[0][2]=0
  copy[0][3]=0
  copy[1][0]=0
  copy[2][0]=0
  copy[3][0]=0
  return Qobj(copy, dims=[[2, 2], [2, 2]])

def measurement(densitymat, c):
  v = 1/(1+c)
  w = 1/(1+2*c)
  x = (c**2+c+1)/(c**2+2*c+1)
  y = (c**2+1)/(c**2+2*c+1)
  z = (2*c**2+np.sqrt(2)*c+1)/(2*c**2+3*c+1)
  copy = np.array(densitymat)
  copy[0][0]= copy[0][0]
  copy[0][1]= copy[0][1]*v
  copy[0][2]= copy[0][2]*v
  copy[0][3]= copy[0][3]*w

  copy[1][0]= copy[1][0]*v
  copy[1][1]= copy[1][1]*x
  copy[1][2]= copy[1][2]*y
  copy[1][3]= copy[1][3]*z

  copy[2][0]= copy[2][0]*v
  copy[2][1]= copy[2][1]*y
  copy[2][2]= copy[2][2]*x
  copy[2][3]= copy[2][3]*z

  copy[3][0]= copy[3][0]*w
  copy[3][1]= copy[3][1]*z
  copy[3][2]= copy[3][2]*z
  copy[3][3]= copy[3][3]
  return Qobj(copy, dims=[[2, 2], [2, 2]])

def fidelity_test(densitymat, ideal_state):
  if ideal_state.type == "ket":
    ideal_state_ket = ideal_state
    ideal_state_bra = ideal_state.dag()
  else:
    ideal_state_bra = ideal_state
    ideal_state_ket = ideal_state.dag()
  # print(ideal_state_bra)
  # print(densitymat)
  # print(ideal_state_ket)
  converted = ideal_state_bra*densitymat*ideal_state_ket
  print(converted)
  return converted[0][0]

initial_state = a0
initial_state = a1*a1.dag()
rho_orig = tensor(initial_state, initial_state)

new_rho = rho_orig

# for i in range(10):
#   op = tensor(r_theta2(np.pi/4, np.pi/10), r_theta2(np.pi/4, np.pi/10))
#   new_rho = unitary_operator(op, new_rho)
#   new_rho = ideal_measurement(new_rho)
#   fig, ax = hinton(new_rho)

ideal_state = Qobj([[0],[np.sqrt(2)/2],[-np.sqrt(2)/2],[0]])
densitymat = ideal_state*ideal_state.dag()
densitymat = Qobj([[0.1,0.5,0.5,0.5],[0.4,0.35,-0.35,0.5],[0.4,-0.35,0.35,0.3],[0.4,0.5,-0.5,0.4]])
# print(ideal_state)
# print(densitymat)
# print('test')
# print(Qobj(fidelity_test(densitymat, ideal_state)))

ideal_state = Qobj([[0],[np.sqrt(2)/2],[-np.sqrt(2)/2],[0]])
c_axis = np.array([])
angle_axis = np.array([])
fidelity_axis = np.array([])
# print(ideal_state)
current_max_fidelity = 0
ideal_angle = [0,0,0] # [c, best angle, fidelity] 
for c in range(1000): #cooperativity
  current_best_angle =[0,0] #[angle, fidelity]
  for i in range(100): #rotation angle
    rotation_angle = np.pi/(i+1)
    new_rho = rho_orig
    current_max_fidelity = 0
    for j in range(i+1): #iterates through rotations
      op = tensor(r_theta2(np.pi, rotation_angle), r_theta2(0, rotation_angle))
      new_rho = unitary_operator(op, new_rho)
      new_rho = Qobj(measurement(new_rho, c+1))
      temp = Qobj(np.array(new_rho), dims=[[4], [4]])
      # print(temp)
      # print(fidelity(ideal_state, temp))
      if fidelity(temp, ideal_state) > current_max_fidelity:
        current_max_fidelity = fidelity(Qobj(np.array(new_rho), dims=[[4], [4]]), ideal_state)
    if current_max_fidelity > current_best_angle[1]: 
     current_best_angle = [rotation_angle, current_max_fidelity]
  ideal_angle = [c+1, current_best_angle[0], current_max_fidelity]
  print(ideal_angle)
  c_axis = np.append(c_axis, ideal_angle[0])
  angle_axis = np.append(angle_axis, ideal_angle[2])
  fidelity_axis = np.append(fidelity_axis, ideal_angle[2])
  # colors = np.append(colors, ideal_angle[2])
# print(x)
# print(y)
x = c_axis
y = fidelity_axis

file_x = open("x_axis.txt","w")
file_y = open("y_axis.txt","w")

for element in x:
    file_x.write(str(element) + "\n")
file_x.close()

for element in y:
    file_y.write(str(element) + "\n")
file_y.close()

plt.subplot(2,1,2)
plt.xscale('log')

a, b, c, d = np.polyfit(x, y, 3)
plt.plot(x, a*x**3+b*x**2+c*x+d)


# plt.scatter(x, y, c=colors, cmap='winter')
# plt.colorbar()'
# plt.plot(x, 1-(1/(x)**(1/3)))
plt.scatter(x, y, color = '#88c999')
plt.xlabel("Cooperativity")
plt.ylabel("Fidelity")

print(a)
print(b)
print(c)
print(d)

plt.show()