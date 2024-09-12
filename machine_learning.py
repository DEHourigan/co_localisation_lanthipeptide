import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, roc_auc_score
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.utils import resample
import numpy as np
import hashlib

# Read in the dataframe 
df = pd.read_csv('/data/san/data0/users/david/intelligence/tables/ctrl_v_LanM_vs_Pks.tsv', sep='\t')
print("Columns in the DataFrame:", df.columns)

# Convert 'NA' to None for proper handling
df['pfam'] = df['pfam'].replace('NA', None)

# Create a pivot table with 1-0 for each pfam, this matrix is used for 1 hot encoding
pivot_df = df.pivot_table(index='Nucleotide_acc', columns='pfam', aggfunc='size', fill_value=0)
pivot_df.reset_index(inplace=True)

# Merge with the original dataframe to get the 'group' labels back
df_group = df[['Nucleotide_acc', 'group']].drop_duplicates()
merged_df = pd.merge(pivot_df, df_group, on='Nucleotide_acc')

# Encoding the 'group' column
le = LabelEncoder()
merged_df['group'] = le.fit_transform(merged_df['group'])

# Define features and labels
X = merged_df.drop(columns=['Nucleotide_acc', 'group'])
y = merged_df['group']

# Function to hash a row
def hash_row(row):
    return hashlib.md5(row.values.tobytes()).hexdigest()

# Hash each row and remove duplicates based on the hashes
X['hash'] = X.apply(hash_row, axis=1)
merged_df['hash'] = X['hash']

# Remove duplicates
merged_df_filtered = merged_df.drop_duplicates(subset=['hash'])

# Redefine features and labels after filtering
X_filtered = merged_df_filtered.drop(columns=['Nucleotide_acc', 'group', 'hash'])
y_filtered = merged_df_filtered['group']

# Identify the larger group
group_counts = y_filtered.value_counts()
larger_group = group_counts.idxmax()
smaller_group_size = group_counts.min()

# Subsample from the larger group to balance the dataset
X_larger_group = X_filtered[y_filtered == larger_group]
y_larger_group = y_filtered[y_filtered == larger_group]
X_smaller_group = X_filtered[y_filtered != larger_group]
y_smaller_group = y_filtered[y_filtered != larger_group]

X_larger_group_subsampled, y_larger_group_subsampled = resample(X_larger_group, y_larger_group, n_samples=smaller_group_size, random_state=42)

# Combine the subsampled larger group with the smaller group
X_balanced = pd.concat([X_larger_group_subsampled, X_smaller_group])
y_balanced = pd.concat([y_larger_group_subsampled, y_smaller_group])

# Train a Random Forest classifier to determine feature importances
clf_initial = RandomForestClassifier(n_estimators=100, random_state=42)
clf_initial.fit(X_balanced, y_balanced)

# Determine feature importance
initial_feature_importances = clf_initial.feature_importances_
initial_feature_names = X_balanced.columns

# Create a DataFrame for the feature importances
initial_importance_df = pd.DataFrame({
    'Feature': initial_feature_names,
    'Importance': initial_feature_importances
})

# Sort the DataFrame by importance
initial_importance_df = initial_importance_df.sort_values(by='Importance', ascending=False)

# Identify the top 50 most important features
top_50_features = initial_importance_df.head(50)['Feature'].tolist()

# Remove the top 50 most important features from the dataset
X_balanced_filtered = X_balanced.drop(columns=top_50_features)

# Split the filtered data into training, validation, and testing sets
X_train_filtered, X_temp_filtered, y_train, y_temp = train_test_split(X_balanced_filtered, y_balanced, test_size=0.3, random_state=42)
X_val_filtered, X_test_filtered, y_val, y_test = train_test_split(X_temp_filtered, y_temp, test_size=0.5, random_state=42)

# Train a Random Forest classifier on the filtered dataset
clf_filtered = RandomForestClassifier(n_estimators=100, random_state=42)
clf_filtered.fit(X_train_filtered, y_train)

# Make predictions on the test data
y_pred_filtered_test = clf_filtered.predict(X_test_filtered)
y_pred_filtered_val = clf_filtered.predict(X_val_filtered)

# Evaluate the model
accuracy_filtered_test = accuracy_score(y_test, y_pred_filtered_test)
accuracy_filtered_val = accuracy_score(y_val, y_pred_filtered_val)
report_filtered_test = classification_report(y_test, y_pred_filtered_test)
report_filtered_val = classification_report(y_val, y_pred_filtered_val)

print(f"Accuracy on test set after filtering: {accuracy_filtered_test}")
print("Classification Report on test set after filtering:")
print(report_filtered_test)

print(f"Accuracy on validation set after filtering: {accuracy_filtered_val}")
print("Classification Report on validation set after filtering:")
print(report_filtered_val)

# Determine feature importance on the filtered dataset
filtered_feature_importances = clf_filtered.feature_importances_
filtered_feature_names = X_balanced_filtered.columns

# Create a DataFrame for the feature importances on the filtered dataset
filtered_importance_df = pd.DataFrame({
    'Feature': filtered_feature_names,
    'Importance': filtered_feature_importances
})

# Sort the DataFrame by importance
filtered_importance_df = filtered_importance_df.sort_values(by='Importance', ascending=False)

# Save the feature importance table
filtered_importance_df.to_csv('filtered_feature_importance.csv', index=False)

# Print the most important features after filtering
print("Top 10 most important features after filtering:")
print(filtered_importance_df.head(10))

# Confusion matrix
cm = confusion_matrix(y_test, y_pred_filtered_test)
plt.figure(figsize=(10, 7))
sns.heatmap(cm, annot=True, fmt='g')
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Confusion Matrix - Test Set')
plt.savefig('confusion_matrix_test.png', dpi=600)
plt.show()

# Bootstrap replicates for the whole set
n_iterations = 10
n_size = int(len(X_balanced_filtered) * 0.7)
bootstrap_accuracies = []

for _ in range(n_iterations):
    # Create a bootstrap sample
    X_sample, y_sample = resample(X_balanced_filtered, y_balanced, n_samples=n_size, random_state=42)
    # Train the model
    clf_bootstrap = RandomForestClassifier(n_estimators=100, random_state=42)
    clf_bootstrap.fit(X_sample, y_sample)
    # Evaluate the model on the test set
    y_pred_bootstrap = clf_bootstrap.predict(X_test_filtered)
    accuracy = accuracy_score(y_test, y_pred_bootstrap)
    bootstrap_accuracies.append(accuracy)

# Plot accuracy distribution
plt.figure(figsize=(10, 7))
plt.hist(bootstrap_accuracies, bins=30, edgecolor='k', alpha=0.7)
plt.axvline(np.mean(bootstrap_accuracies), color='r', linestyle='dashed', linewidth=1)
plt.title('Bootstrap Accuracy Distribution')
plt.xlabel('Accuracy')
plt.ylabel('Frequency')
plt.savefig('bootstrap_accuracy_distribution.png')
plt.show()

# ROC Curve
y_pred_prob = clf_filtered.predict_proba(X_test_filtered)[:, 1]
fpr, tpr, _ = roc_curve(y_test, y_pred_prob)
roc_auc = roc_auc_score(y_test, y_pred_prob)

plt.figure(figsize=(5, 5))
plt.plot(fpr, tpr, color='red', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='gold', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")
plt.savefig('roc_curve.png', dpi=600)
plt.show()
