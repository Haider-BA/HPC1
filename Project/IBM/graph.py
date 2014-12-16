import matplotlib.pyplot as plt
import numpy as np

#make data
u_x = np.array([-1, 0, 1])
u_y = np.array([0., 0., 0.])

v_x = np.array([-0.5, 0.5, -0.5, 0.5])
v_y = np.array([-0.5, -0.5, 0.5, 0.5,])

p_x = np.array([-0.5, 0.5])
p_y = np.array([0, 0])

#plot marks
plt.plot(u_x, u_y, marker=r'$\bullet $', markersize=8, linewidth=0, alpha=0.6, color='b', label='u_flux')

plt.plot(v_x, v_y, marker=r'$\bullet $', markersize=8, linewidth=0, alpha=0.5, color='g', label='v_flux')

plt.plot(p_x, p_y, marker=r'$\circ $', markersize=8, linewidth=0, color ='r', label='p_flux')

#plot annotation
plt.annotate('$i-1$', xy=(-1.5, -1.0), xytext=(+5,+5), textcoords='offset points')
plt.annotate('$i$', xy=(-0.5, -1.0), xytext=(+5,+5), textcoords='offset points')
plt.annotate('$i+1$', xy=(0.5, -1.0), xytext=(+5,+5), textcoords='offset points')
plt.annotate('$i+2$', xy=(1.5, -1.0), xytext=(+5,+5), textcoords='offset points')

plt.annotate('$j-1$', xy=(-1.5, -1.0), xytext=(+5,+5), textcoords='offset points')
plt.annotate('$j$', xy=(-1.5, 0), xytext=(+5,+5), textcoords='offset points')
plt.annotate('$j+1$', xy=(-1.5, 1.0), xytext=(+5,+5), textcoords='offset points')

plt.annotate('$u_{i-1,j}$', xy=(-1, 0), xytext=(+5,+5), textcoords='offset points')
plt.annotate('$u_{i,j}$', xy=(0, 0), xytext=(+5,+5), textcoords='offset points')
plt.annotate('$u_{i+1,j}$', xy=(1,0), xytext=(+5,+5), textcoords='offset points')

plt.annotate('$v_{i,j-1}$', xy=(-0.5, -0.5), xytext=(+5,+5), textcoords='offset points')
plt.annotate('$v_{i+1,j-1}$', xy=(0.5, -0.5), xytext=(+5,+5), textcoords='offset points')
plt.annotate('$v_{i,j}$', xy=(-0.5, 0.5), xytext=(+5,+5), textcoords='offset points')
plt.annotate('$v_{i+1,j}$', xy=(0.5, 0.5), xytext=(+5,+5), textcoords='offset points')

plt.annotate('$p_{i,j}$', xy=(-0.5, 0), xytext=(+5,+5), textcoords='offset points')
plt.annotate('$p_{i+1,j}$', xy=(0.5, 0), xytext=(+5,+5), textcoords='offset points')

plt.grid(True, color='k', linewidth=2)

#plot shadow
shadow_u_x = np.array([-0.5, 0.5, 0.5, -0.5])
shadow_u_y = np.array([-0.5, -0.5, 0.5, 0.5])
shadow_v_x = np.array([-1.0, 0.0, 0.0, -1.0])
shadow_v_y = np.array([0.0, 0.0,  1.0, 1.0])
shadow_p_x = np.array([-1.0, 0.0, 0.0, -1.0])
shadow_p_y = np.array([-0.5, -0.5, 0.5, 0.5])

plt.fill(shadow_u_x, shadow_u_y, linewidth=0, alpha=0.2, color='b')
plt.fill(shadow_v_x, shadow_v_y, linewidth=0, alpha=0.2, color='g')
plt.fill(shadow_p_x, shadow_p_y, linewidth=0, alpha=0.2, color='r')

#set the figure style
ax = plt.gca()
#ax.set_xticks(np.arange(-2,2,0.5))
#ax.set_yticks(np.arange(-2,2,0.5))

ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['top'].set_color('none')

plt.axis('scaled')
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)

#add the legend
ax.legend(loc='right center')

plt.show()



