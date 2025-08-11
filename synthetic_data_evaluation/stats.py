# import matplotlib.pyplot as plt

def printStats(depths):

    # After the loop
    n = len(depths)
    maxD = max(depths)
    mean = sum(depths) / n
    variance = sum((d - mean) ** 2 for d in depths) / (n - 1)  # sample stddev
    stddev = variance ** 0.5

    # plt.hist(depths, bins=100, edgecolor='black')
    # plt.xlabel('Count')
    # plt.ylabel('Frequency')
    # plt.title('Histogram of Counts')
    # plt.grid(True)
    # plt.show()

    print(f"Average: {float(mean):.3f}")
    print(f"Max: {float(maxD):.3f}")
    print(f"Standard deviation: {float(stddev):.3f}")
