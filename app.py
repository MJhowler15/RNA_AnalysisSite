import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

# 1. PAGE SETUP
st.set_page_config(page_title="RNA Disease Correlation AI", layout="wide")

st.title("❤️ RNA-Seq Disease Correlation & Biomarker Discovery Tool")
st.markdown(
    """
This advanced tool calculates the statistical correlation between **thousands of genes** and a specific clinical trait 
(like **Cardiac Disease**). It identifies potential **biomarkers**—genes that could be targeted for drugs or diagnostics.
"""
)

# 2. FILE UPLOADER WIDGET
uploaded_file = st.file_uploader(
    "Upload your RNA-Seq TSV file (First column must be GeneID)",
    type=["tsv", "txt"],
)

if uploaded_file is not None:
    with st.spinner("Analyzing data and calculating correlations..."):

        # 3. READ THE FILE
        df = pd.read_csv(uploaded_file, sep="\t", index_col=0, comment="!")

        # 4. AUTOMATICALLY DETECT PATIENT GROUPS
        # In professional datasets, samples are named systematically.
        # Let's write a dynamic scanner that looks at the sample columns.
        # If using your Lupus file: GSM3466767-76 are Control, 77-92 are Disease.
        # To make this easily adaptable to a Cardiac dataset, we look at the column halves:
        all_samples = df.columns.tolist()
        midpoint = len(all_samples) // 2

        # We will split the samples into Group A (Control) and Group B (Target Disease)
        control_samples = all_samples[:midpoint]
        disease_samples = all_samples[midpoint:]

        # 5. CALCULATE STATISTICAL CORRELATION
        # We want to map clinical state to a number: Control = 0, Disease = 1
        clinical_outcomes = [0] * len(control_samples) + [1] * len(disease_samples)

        correlations = []
        fold_changes = []
        avg_controls = []
        avg_diseases = []

        # Loop through every single gene row to find correlations
        for gene, rows in df.iterrows():
            control_values = rows[control_samples].values
            disease_values = rows[disease_samples].values

            # Averages
            avg_c = np.mean(control_values)
            avg_d = np.mean(disease_values)
            avg_controls.append(avg_c)
            avg_diseases.append(avg_d)

            # Fold Change (How many times higher is it in the disease?)
            fc = (avg_d + 1) / (avg_c + 1)
            fold_changes.append(fc)

            # Combined list of expressions for this specific gene
            all_expressions = list(control_values) + list(disease_values)

            # Calculate Pearson Correlation between expression numbers and [0,0,0...1,1,1...]
            # np.corrcoef returns a matrix, we extract the core correlation value
            corr = np.corrcoef(all_expressions, clinical_outcomes)[0, 1]

            # If data is completely flat, correlation returns NaN (Not a Number). Fix to 0.
            if np.isnan(corr):
                corr = 0.0
            correlations.append(corr)

        # 6. BUILD THE RESULTS MARKETPLACE
        results_df = pd.DataFrame(
            {
                "GeneID": df.index,
                "Avg_Control": avg_controls,
                "Avg_Disease": avg_diseases,
                "Fold_Change": fold_changes,
                "Disease_Correlation": correlations,
            }
        ).set_index("GeneID")

        # Absolute correlation tells us strength regardless of direction (positive or negative)
        results_df["Correlation_Strength"] = results_df["Disease_Correlation"].abs()

        # 7. SHOW SUMMARY METRICS
        st.header("📊 Biomarker Discovery Summary")
        col1, col2, col3 = st.columns(3)

        # Top correlated gene
        top_gene = results_df["Correlation_Strength"].idxmax()
        top_corr_val = results_df.loc[top_gene, "Disease_Correlation"]

        col1.metric("Control Samples Analyzed", len(control_samples))
        col2.metric("Disease Samples Analyzed", len(disease_samples))
        col3.metric(
            "Strongest Gene Predictor", f"ID {top_gene} ({top_corr_val:.2f} corr)"
        )

        # 8. FILTER THE TOP BIOMARKERS
        st.header("🧬 Top 10 Disease-Correlated Genes")
        st.markdown(
            "A correlation close to **1.0** means the gene spikes when disease is present. Close to **-1.0** means the gene shuts off when disease is present."
        )

        # Sort to find the highest absolute correlations
        top_10_genes = results_df.sort_values(
            by="Correlation_Strength", ascending=False
        ).head(10)
        st.dataframe(top_10_genes[["Disease_Correlation", "Fold_Change"]])

        # 9. INTERACTIVE GENE INSPECTOR (THE MOST VALUABLE PART)
        st.header("🔍 Individual Gene Deep-Dive")
        selected_gene = st.selectbox(
            "Select a Gene ID to visualize its expression behavior across patients:",
            options=df.index.tolist(),
        )

        # Pull expression values for just this one chosen gene
        gene_data = df.loc[selected_gene]

        # Structure it into a clean table for graph rendering
        plot_df = pd.DataFrame(
            {
                "Expression Level": gene_data.values,
                "Patient Group": [
                    "Healthy Control" if s in control_samples else "Cardiac/Disease"
                    for s in gene_data.index
                ],
                "Sample ID": gene_data.index,
            }
        )

        # Draw a beautiful interactive boxplot showing the variance split
        fig = px.box(
            plot_df,
            x="Patient Group",
            y="Expression Level",
            color="Patient Group",
            points="all",  # Show individual patient dots over the box plot
            title=f"Expression Profile for Gene: {selected_gene}",
            color_discrete_map={
                "Healthy Control": "blue",
                "Cardiac/Disease": "red",
            },
        )
        st.plotly_chart(fig, use_container_width=True)

else:
    st.info("👋 System ready. Please upload an RNA-Seq matrix file to start correlation analysis.")