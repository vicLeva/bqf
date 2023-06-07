import matplotlib.pyplot as plt

x = [-0.875 + 2*z for z in range(0,9)]
print(x)
plt.plot(x, 'bo')
plt.show()