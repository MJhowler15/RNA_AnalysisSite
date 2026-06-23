import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# 1. SETUP THE WEB PAGE INTERFACE
st.set_page_config(page_title="AI RNA-Seq Explorer", layout="wide")

st.title("🤖 Adaptive AI RNA-Seq Data Explorer")
st.markdown(
    """
Upload **any** RNA-seq expression matrix file (TSV format). 
The AI will automatically handle the preprocessing, perform **Principal Component Analysis (PCA)** to find hidden patterns, 
and use **K-Means Machine Learning** to cluster similar samples together.
"""
)

# 2. CREATE A DRAG-AND-DROP FILE UPLOADER WIDGET
uploaded_file = st.file_uploader(
    "Choose an RNA-Seq TSV file (The first column must be GeneID)", type=["tsv"]
)

# If the user has uploaded a file, start the AI magic
if uploaded_file is not None:
    # Tell the user the app is working
    with st.spinner("AI is processing your dataset..."):

        # 3. READ THE UPLOADED FILE
        # We read the file, setting the 'GeneID' column as our index rows.
        df = pd.read_csv(uploaded_file, sep="\t", index_col=0)

        # 4. TRANSPOSE DATA FOR MACHINE LEARNING
        # Machine learning expects samples as rows and genes as features.
        # Original shape: Rows = 20,000+ genes, Columns = 26 samples.
        # Transposed shape: Rows = 26 samples, Columns = 20,000+ genes.
        df_transposed = df.T

        # 5. NORMALIZE AND SCALE DATA
        # RNA-seq data can range from 0 to millions. ML models break if the scale is too wide.
        # StandardScaler shifts data so the mean is 0 and variance is 1.
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(df_transposed)

        # 6. APPLICATION OF MACHINE LEARNING (PCA)
        # We compress the 20,000+ gene columns down to just 2 'Principal Components' (patterns)
        pca = PCA(n_components=2)
        pca_results = pca.fit_transform(scaled_data)

        # Create a neat dataframe of our AI dimensions
        pca_df = pd.DataFrame(
            data=pca_results,
            columns=["PC1 (Main Pattern)", "PC2 (Secondary Pattern)"],
        )
        # Keep track of the actual Sample IDs (like GSM3466767)
        pca_df["SampleID"] = df_transposed.index

        # 7. APPLICATION OF MACHINE LEARNING (K-MEANS CLUSTERING)
        st.sidebar.header("AI Model Controls")
        # Let the user choose how many distinct groups the AI should look for
        num_clusters = st.sidebar.slider(
            "How many biological groups do you expect?",
            min_value=2,
            max_value=5,
            value=2,
        )

        # Initialize the K-Means AI model
        kmeans = KMeans(n_clusters=num_clusters, random_state=42)
        # Train the AI model to group the sample rows
        pca_df["AI_Cluster_Group"] = kmeans.fit_predict(scaled_data)
        # Convert group numbers (0, 1) to text labels for the graph legend
        pca_df["AI_Cluster_Group"] = pca_df["AI_Cluster_Group"].astype(str)

        # 8. PRESENTING RESULTS: METRICS SECTION
        st.header("📊 Dataset Overview")
        col1, col2, col3 = st.columns(3)
        col1.metric("Total Samples Found", len(df_transposed))
        col2.metric("Total Genes Analyzed", len(df))
        # Information retention metric tells us how much biological truth our 2D compression kept
        variance_retained = (pca.explained_variance_ratio_.sum()) * 100
        col3.metric("AI Info Retention", f"{variance_retained:.1f}%")

        # 9. PRESENTING RESULTS: THE INTERACTIVE GRAPH
        st.header("📈 AI Unsupervised Sample Clustering")
        st.markdown(
            "Each dot is a human sample. Dots closer together have highly similar gene expression behaviors."
        )

        fig = px.scatter(
            pca_df,
            x="PC1 (Main Pattern)",
            y="PC2 (Secondary Pattern)",
            color="AI_Cluster_Group",  # Color dots automatically by what the AI discovered
            text="SampleID",  # Put the text name next to the dot
            title="PCA & K-Means Mapping",
            labels={"AI_Cluster_Group": "AI-Discovered Group"},
        )
        fig.update_traces(marker=dict(size=12))
        st.plotly_chart(fig, use_container_width=True)

        # 10. PRESENTING RESULTS: RAW DATA LOOKUP
        st.header("📋 Raw Expression Sub-Matrix")
        st.write("First 50 genes across all your samples:")
        st.dataframe(df.head(50), use_container_width=True)

else:
    # If the user hasn't uploaded a file yet, show a welcoming notice
    st.info("👋 Welcome! Please upload your RNA-seq dataset (.tsv) file to begin.")