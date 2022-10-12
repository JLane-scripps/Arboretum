import streamlit as st

st.title('Arborist Monitor')

arborist_save_files = st.file_uploader("PSM pickle files", type=['.pkl'], accept_multiple_files=True)

col1, col2 = st.columns(2)
mz = col1.number_input(label='mz', min_value=0.0, step=0.01)
ppm = col2.number_input(label='ppm', min_value=0.0, step=0.01)

col1, col2 = st.columns(2)
rt = col1.number_input(label='rt', min_value=0.0, step=0.01)
offset = col2.number_input(label='offset', min_value=0.0, step=0.01)

col1, col2 = st.columns(2)
ook0 = col1.number_input(label='ook0', min_value=0.0, step=0.01)
tolerance = col2.number_input(label='tolerance', min_value=0.0, step=0.01)

# Add other input

if st.button("Run"):
    for file in arborist_save_files:
        print(file.name)
        # load into PSM tree


