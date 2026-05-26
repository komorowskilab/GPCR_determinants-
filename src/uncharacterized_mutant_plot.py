import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# 1. Data
# -----------------------------
data1 = {
    "Sample": ["M21","M22","M23","M26","M27", "M42",
                "M48","M47","M50","M51","M53","M56",
                "M64","M78","M79","M83","M84","M85",
                "M87","M89"],
    "Adrenergic": [0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05],
    "Octopamine": [0.04494382,
0.04494382,
0.04494382,
0.05617978,
0.04494382,
0.06741573,
0.06741573,
0.06741573,
0.04494382,
0.05617978,
0.04494382,
0.05617978,
0.05617978,
0.04494382,
0.04494382,
0.04494382,
0.04494382,
0.04494382,
0.04494382,
0.03370787]
}



data2 = {
    "Sample": ["M21","M22","M23","M26","M27", "M42",
                "M48","M47","M50","M51","M53","M56",
                "M64","M78","M79","M83","M84","M85",
                "M87","M89"],
    "Adrenergic": [0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05,
0.05], 
"Octopamine": [0.05172414,
0.05172414,
0.03448276,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.05172414,
0.03448276]
}


data3 = {
    "Sample": ["M21","M22","M23","M26","M27", "M42",
                "M48","M47","M50","M51","M53","M56",
                "M64","M78","M79","M83","M84","M85",
                "M87","M89"],
    "Adrenergic": [
0.012578616,
0.006289308,
0.081761006,
0.062893082,
0.037735849,
0.037735849,
0.037735849,
0.075471698,
0.050314465,
0.069182390,
0.044025157,
0.069182390,
0.056603774,
0.044025157,
0.075471698,
0.075471698,
0.075471698,
0.018867925,
0.037735849,
0.031446541],
"Octopamine": [
0.07042254,
0.08450704,
0.02816901,
0.11267606,
0.07042254,
0.04225352,
0.02816901,
0.02816901,
0.04225352,
0.01408451,
0.01408451,
0.07042254,
0.04225352,
0.04225352,
0.02816901,
0.02816901,
0.05633803,
0.08450704,
0.05633803,
0.05633803]
}



data4 = {
    "Sample": ["M21","M22","M23","M26","M27", "M42",
                "M48","M47","M50","M51","M53","M56",
                "M64","M78","M79","M83","M84","M85",
                "M87","M89"],
    "Adrenergic": [
0.02010050,
0.01507538,
0.07537688,
0.06030151,
0.04020101,
0.04020101,
0.04020101,
0.07035176,
0.05025126,
0.06532663,
0.04522613,
0.06532663,
0.05527638,
0.04522613,
0.07035176,
0.07035176,
0.07035176,
0.02512563,
0.04020101,
0.03517588],
"Octopamine": [
0.05911330,
0.06403941,
0.02463054,
0.07389163,
0.05911330,
0.04926108,
0.04433498,
0.04433498,
0.04926108,
0.03940887,
0.03940887,
0.05911330,
0.04926108,
0.04926108,
0.04433498,
0.04433498,
0.05418719,
0.06403941,
0.05418719,
0.03448276]
}



df1 = pd.DataFrame(data1).set_index("Sample")
df2 = pd.DataFrame(data2).set_index("Sample")
df3 = pd.DataFrame(data3).set_index("Sample")
df4 = pd.DataFrame(data4).set_index("Sample")

models = [df1, df2, df3, df4]
cols = ["Adrenergic", "Octopamine"]

# -----------------------------
# 2. Normalize each model (row-wise)
# -----------------------------
normalized_models = []

for df in models:
    row_sums = df[cols].sum(axis=1)
    df_norm = df.copy()
    df_norm[cols] = df[cols].div(row_sums, axis=0)
    normalized_models.append(df_norm)

# -----------------------------
# 3. Combine models (average)
# -----------------------------
combined = sum(df[cols] for df in normalized_models) / len(normalized_models)

# -----------------------------
# 4. Predicted class
# -----------------------------
combined["predictedClass"] = combined.idxmax(axis=1)

print("Combined normalized votes:")
print(combined)

# -----------------------------
# 5. Plot
# -----------------------------
combined[cols].plot(kind="bar", stacked=True)

plt.ylabel("Normalized vote")
plt.title("Votes per sample (combined models)")
plt.xticks(rotation=0)
plt.tight_layout()
plt.show()

# -----------------------------
# 5. Plot with 0.5 bar highlighted
# -----------------------------
ax = combined[cols].plot(kind="bar", stacked=True)

plt.ylabel("Normalized vote")
plt.title("Votes per sample (combined models)")
plt.xticks(rotation=0)

# Option 1: Add a vertical line at 0.5
plt.axhline(y=0.5, color='red', linestyle='--', linewidth=2, label='0.5 threshold')

# Option 2: Annotate the 0.5 line
plt.text(-0.5, 0.5, '0.5', fontsize=12, color='red', 
         verticalalignment='center', fontweight='bold')

plt.legend()
plt.tight_layout()
plt.show()