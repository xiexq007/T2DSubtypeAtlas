
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, jaccard_score
from sklearn.utils import resample

# --- Configuration & Constants ---

FILE_PATH = "./43_donorDataForCluster.csv"
CLUSTER_FEATURES = ['bmi', 'age', 'hba1c', 'c_peptide'] 
N_CLUSTERS = 4
RANDOM_STATE = 42
COLORS = ['#367eb6', '#4aaf49', '#e31a1a', '#fea502']

# Mapping dictionary for cluster labeling based on results
FEMALE_MAPPING = {0: 'MOD', 1: 'MARD', 2: 'SIDD', 3: 'SIRD'}
MALE_MAPPING = {0: 'MOD', 1: 'SIDD', 2: 'MARD', 3: 'SIRD'}

def perform_clustering(df, features, n_clusters, random_state=42):
   
    """
    Scales features and performs KMeans clustering on the provided dataframe.
    """
    
    scaler = MinMaxScaler()
    X = df[features]
    X_scaled = scaler.fit_transform(X)
    
    kmeans = KMeans(
        n_clusters=n_clusters, 
        init='k-means++', 
        n_init=1000, 
        random_state=random_state
    )
    labels = kmeans.fit_predict(X_scaled)
    score = silhouette_score(X_scaled, labels)
    
    return labels, score

def plot_cluster_feature_box(data, x, y, colors, ylabel, figsize=(4, 4)):
    
    """
    Generates a stylized boxplot for clinical features across clusters.
    """

    target_order = ['MARD', 'MOD', 'SIDD', 'SIRD']
    data[x] = pd.Categorical(data[x], categories=target_order, ordered=True)
    
    plt.figure(figsize=figsize)
    ax = sns.boxplot(
        x=x, y=y, data=data, hue=x, 
        palette=colors, width=0.65, linewidth=1.8, legend=False
    )
    
    plt.ylabel(ylabel)
    plt.xlabel('')
    
    # Customize spine thickness
    for spine in ax.spines.values():
        spine.set_linewidth(2)
        
    plt.tight_layout()
    plt.show()

def main():
    
    # 1. Load data
    data = pd.read_csv(FILE_PATH)
    
    # 2. Split data by gender
    data_female = data[data['gender'] == 'Female'].copy()
    data_male = data[data['gender'] == 'Male'].copy()

    # 3. Perform Clustering for each group
    features = data.columns[2:6].tolist() 
    
    f_labels, f_score = perform_clustering(data_female, features, N_CLUSTERS, RANDOM_STATE)
    m_labels, m_score = perform_clustering(data_male, features, N_CLUSTERS, RANDOM_STATE)
    
    print(f"Female Silhouette Score: {f_score:.3f}")
    print(f"Male Silhouette Score: {m_score:.3f}")

    # 4. Assign labels and map to clinical subtypes
    data_female['kclusterGroup'] = pd.Series(f_labels, index=data_female.index).map(FEMALE_MAPPING)
    data_male['kclusterGroup'] = pd.Series(m_labels, index=data_male.index).map(MALE_MAPPING)

    # 5. Recombine data
    data_combined = pd.concat([data_female, data_male], ignore_index=True)

    # 6. Visualization (Figure 1A-D) ----------
    plot_configs = [
        ('age', 'Age (years)'),
        ('bmi', 'BMI (kg/m²)'),
        ('c_peptide', 'C-peptide (ng/ml)'),
        ('hba1c', 'HbA1c (%)'),
    ]

    for column, label in plot_configs:
        plot_cluster_feature_box(
            data=data_combined, 
            x='kclusterGroup', 
            y=column, 
            colors=COLORS, 
            ylabel=label
        )
    
   

if __name__ == "__main__":
    main()

# Figure 1E ----------

# Data extracted from the clustering results
values = [17, 13, 10, 3]  
labels = ['MARD', 'MOD', 'SIDD', 'SIRD']  
colors = ['#4aaf49','#367eb6','#fea502','#e31a1a'] 

plt.figure(figsize=(5, 5))
plt.pie(
    values, 
    labels=labels, 
    colors=colors, 
    autopct='%1.1f%%', 
    startangle=140)

plt.axis('equal')  
plt.show()

# Figure 1F ----------
# Data Prearation Functions

def get_normalized_gender_df(data_combined, gender):
    
    """
    Filters data by gender, normalizes clinical features, 
    and returns a structured DataFrame with donor_ID.
    """
    
    # Select columns: donor_ID (0), gender (1), and features (2:6)
    subset = data_combined.iloc[:, [0, 1, 2, 3, 4, 5]]
    filtered = subset[subset['gender'] == gender].drop(columns=['gender'])
    
    scaler = MinMaxScaler()
    
    features_to_scale = filtered.iloc[:, 1:]
    scaled_data = scaler.fit_transform(features_to_scale)
    
    normalized_df = pd.DataFrame(scaled_data, columns=features_to_scale.columns)
    
    # Reset index to align donor_ID
    final_df = pd.concat([filtered.iloc[:, [0]].reset_index(drop=True), normalized_df], axis=1)
    return final_df

def get_standard_labels(data_combined, gender):
    
    """
    Creates a dictionary mapping Cluster Names to lists of Donor IDs (Ground Truth).
    """
    
    filtered = data_combined[data_combined['gender'] == gender]
    mapping = filtered.groupby('kclusterGroup')['donor_ID'].apply(list).to_dict()
    return mapping

# Robustness Analysis (Jaccard Calculation)
def run_jaccard_bootstrapping(df, df_must_include, standard_labels, feature_to_class, n_iterations=200, n_samples=70):
    
    """
    Performs bootstrap resampling and calculates Jaccard similarity for cluster stability.
    """
   
    results = []
    
    for i in range(n_iterations):
        
        # Bootstrap sampling with replacement (set seeds to marke results fix)
        sampled = resample(df, replace=True, n_samples=n_samples)
        
        # Ensure 'must-include' samples are present (Make sure to randomly select at least one sample for each category before returning to the random sampling)
        combined_sample = pd.concat([sampled, df_must_include], ignore_index=True).drop_duplicates()
        
        # KMeans Clustering
        kmeans = KMeans(n_clusters=4, init='k-means++', n_init=100, random_state=42)
        X = combined_sample.drop(columns=['donor_ID'])
        clusters = kmeans.fit_predict(X)
        
        # Match clusters to clinical subtypes based on feature means
        cluster_means = X.groupby(clusters).mean()
        class_to_donor_ids = {}
        used_clusters = set()

        # Match MOD, SIDD, SIRD based on BMI, HbA1c, C-peptide maximums
        for feature, target_class in feature_to_class.items():
            if feature != 'age':
                best_cluster = cluster_means[feature].idxmax()
                class_to_donor_ids[target_class] = combined_sample[clusters == best_cluster]['donor_ID'].tolist()
                used_clusters.add(best_cluster)

        # Match MARD to the remaining cluster (Age-related)
        remaining = set(range(4)) - used_clusters
        mard_class = feature_to_class['age']
        for c in remaining:
            class_to_donor_ids[mard_class] = combined_sample[clusters == c]['donor_ID'].tolist()

        # Calculate Jaccard Similarity for each group
        iter_scores = {}
        for group in ['MOD', 'MARD', 'SIDD', 'SIRD']:
            set_std = set(standard_labels.get(group, []))
            set_new = set(class_to_donor_ids.get(group, []))
            
            intersection = len(set_std.intersection(set_new))
            union = len(set_std.union(set_new))
            j_index = intersection / union if union != 0 else 0
            iter_scores[group] = j_index
            
        iter_scores['Average'] = np.mean(list(iter_scores.values()))
        results.append(iter_scores)

    return pd.DataFrame(results)

# Visualization
def plot_jaccard_barplot(score_df):
    
    """
    Plots the final Jaccard scores for both genders.
    """
    
    plt.figure(figsize=(8, 5))
    color_dict = {'Male': '#1f78b4', 'Female': '#a6cee3'}
    
    ax = sns.barplot(x='Group', y='Jaccard', hue='Gender', data=score_df, palette=color_dict)
    
    plt.axhline(y=0.75, color='black', linestyle='--', alpha=0.6)
    plt.ylabel('Jaccard Score', fontsize=12)
    plt.xlabel('')
    plt.legend()
    
    # Stylize spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    
    # A. Setup Constants
    FEATURE_MAPPING = {
        'bmi': 'MOD',
        'age_years': 'MARD',
        'hba1c': 'SIDD',
        'c_peptide_ng_ml': 'SIRD'
    }
    MALE_MUST_IDS = ["HPAP-161", "HPAP-108", "HPAP-083", "HPAP-070"]
    FEMALE_MUST_IDS = ["HPAP-013", "HPAP-061", "HPAP-126", "HPAP-145"]

    # B. Processing
    # Get normalized dataframes
    male_df_norm = get_normalized_gender_df(data_combined, 'Male')
    female_df_norm = get_normalized_gender_df(data_combined, 'Female')
    
    # Get standard labels
    male_std_labels = get_standard_labels(data_combined, 'Male')
    female_std_labels = get_standard_labels(data_combined, 'Female')
    
    # Must-include subsets
    male_must_df = male_df_norm[male_df_norm['donor_ID'].isin(MALE_MUST_IDS)]
    female_must_df = female_df_norm[female_df_norm['donor_ID'].isin(FEMALE_MUST_IDS)]

    # C. Run Analysis (Example for Female)
    print("Running stability analysis...")
    female_results = run_jaccard_bootstrapping(
        female_df_norm, female_must_df, female_std_labels, FEATURE_MAPPING
    )
 
    # D. Final Visualization
    jaccard_summary_data = {
        'Group': ['MARD','MARD','MOD','MOD','SIDD','SIDD','SIRD','SIRD','Average','Average'],
        'Gender':['Male', 'Female'] * 5,  
        'Jaccard':[0.846, 0.772, 0.823, 0.750, 0.937, 0.899, 0.771, 0.804, 0.844, 0.806]
    }
    plot_jaccard_barplot(pd.DataFrame(jaccard_summary_data))