import matplotlib.pyplot as plt
data = []
with open('C:\\Users\\NoraS\\Dropbox\\Dokumente\\Studium\\KT\\Praktikum\\Messdaten\\PMT1_energy_nogate - Kopie.Spe') as file:
    for line in file.readlines():
        data.append(line[:-1])

x = range(len(data))
plt.plot(x,data)
print(data[20:40])