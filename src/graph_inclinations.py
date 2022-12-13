import matplotlib.pyplot as plt

def graph():
    plt.hist(incs, bins=tuple(range(145)))
    plt.title('Debris in Low Earth Orbit by Inclination', fontsize=12)
    plt.xlabel('Inclination')
    plt.ylabel('Number of Cataloged Objects')
    plt.yscale('log')
    plt.gcf().set_size_inches(6, 3)
    plt.savefig(r'graphs\leo_inclinations.png', bbox_inches='tight', dpi=300)
    plt.clf()

if __name__ == '__main__':
    with open('../data/incs.txt', mode='r', encoding='UTF-8') as f:
        incs = f.readlines()
    incs = list(map(lambda s: float(s.strip()), incs))
    graph()
