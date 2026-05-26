import matplotlib.pyplot as plt
import numpy as np

score_dict = {}
bad_results = {}

# Read the pLDDT scores
with open("plddt_score_table.tsv", "r") as f:
    for line in f:
        if line.startswith("Sample"):
            continue
        else:
            fields = line.strip().split("\t")
            name = fields[0]
            score = float(fields[1])
            score_dict[name] = score  


# Count samples below threshold
count_below_60 = sum(1 for score in score_dict.values() if score < 60)

# Get bad results (scores < 60)
for name, score in score_dict.items():  
    if score < 60:
        bad_results[name] = score

# Print results
print(f"Total samples: {len(score_dict)}")
print(f"Bad results (pLDDT < 60): {len(bad_results)}")
print("\nSamples with low pLDDT scores:")
for name, score in sorted(bad_results.items(), key=lambda x: x[1]):
    print(f"  {name}: {score:.3f}")

# Create histogram
scores = list(score_dict.values())

plt.figure(figsize=(10, 6))
plt.hist(scores, bins=30, edgecolor='black', alpha=0.7)
plt.axvline(x=60, color='red', linestyle='--', linewidth=2, label='Threshold (60)')
plt.axvline(x=np.mean(scores), color='green', linestyle='--', linewidth=2, label=f'Mean ({np.mean(scores):.2f})')
plt.axvline(x=np.median(scores), color='blue', linestyle='--', linewidth=2, label=f'Median ({np.median(scores):.2f})')

plt.xlabel('pLDDT Score', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.title('Distribution of pLDDT Scores', fontsize=14, fontweight='bold')
plt.legend()
plt.grid(True, alpha=0.3)


# Add statistics text box
stats_text = f'N = {len(scores)}\nMean = {np.mean(scores):.2f}\nMedian = {np.median(scores):.2f}\nStd = {np.std(scores):.2f}\nMin = {min(scores):.2f}\nMax = {max(scores):.2f}\n\nBelow threshold (60) : {count_below_60}'
plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))


plt.tight_layout()
plt.savefig('plddt_distribution.png', dpi=300)
print("\nHistogram saved as 'plddt_distribution.png'")
plt.show()

# Check for skewness
from scipy import stats as scipy_stats
skewness = scipy_stats.skew(scores)
print(f"\nSkewness: {skewness:.3f}")
if abs(skewness) < 0.5:
    print("Distribution is approximately symmetric (normal)")
elif skewness < -0.5:
    print("Distribution is left-skewed (negatively skewed)")
else:
    print("Distribution is right-skewed (positively skewed)")


            