import tkinter as tk
# import customtkinter as ctk 
from tkinter import font
from tkinter import messagebox, scrolledtext
from subprocess import Popen, PIPE, STDOUT
import threading
import subprocess
import os
from tkinter import ttk
import signal

# def run_script():
#     """ Run the script in the same terminal and handle CTRL + C """
#     try:
#         # Define the script to run
#         script_path = os.path.join(os.path.dirname(__file__), "../run.sh")
        
#         #Running the script in the same terminal using Popen
#         process = subprocess.Popen(
#             ["bash", script_path] if os.name == "posix" else ["cmd", "/c", script_path],
#             stdout=subprocess.PIPE,
#             stderr=subprocess.PIPE,
#             text=True
#         )
        
#         # process = subprocess.Popen(
#         #     script_path,
#         #     shell=True,
#         #     stdin=subprocess.PIPE,
#         #     stdout=subprocess.PIPE,
#         #     stderr=subprocess.PIPE,
#         #     universal_newlines=True
#         # )
        
        
#         # Reading the output and error streams
#         for line in process.stdout:
#             print(line, end="")  # Print output in real-time to the terminal
#         for line in process.stderr:
#             print(line, end="")  # Print errors in real-time to the terminal

#         # Wait for the process to complete
#         process.wait()
        

#     except KeyboardInterrupt:
#         print("\nScript interrupted with CTRL + C. You can run the script again.")
#         # Optionally handle cleanup or restart logic here
#     except Exception as e:
#         messagebox.showerror("Error", str(e))
#         print(f"Error: {str(e)}\n")
        
# def run_script():
#     # Path to your bash script
#     script_path = os.path.join(os.path.dirname(__file__), "../run.sh")

#     # Open the bash script with Popen
#     process = subprocess.Popen(
#         script_path,
#         shell=True,
#         stdout=subprocess.PIPE,
#         stderr=subprocess.PIPE,
#         text=True
#     )

#     # Continuously read the output
#     for line in iter(process.stdout.readline, ''):
#         terminal_output.insert(tk.END, line)
#         terminal_output.see(tk.END)  # Scroll to the end of the text widget

#     process.stdout.close()
#     process.wait()

#     # Capture any errors
#     if process.returncode != 0:
#         error_output = process.stderr.read()
#         terminal_output.insert(tk.END, error_output)
#         terminal_output.see(tk.END)
    
#     process.stderr.close() 

    
def run_script_in_terminal():
    # Path to your bash script or Python script
    command = os.path.join(os.path.dirname(__file__), "../run.sh") # or ['python3', 'your_python_script.py']

    # Run the command in the terminal without capturing the output in the GUI
    subprocess.run(command, shell=True)

def start_script():
    # Run the script in a separate thread to keep the GUI responsive
    threading.Thread(target=run_script_in_terminal, daemon=True).start()
    
    

def show_format():
    format_text = (
        "The 'parameters.dat' file should have the following format:\n\n"
        "xlen = <float>\n"
        "ylen = <float>\n"
        "zlen = <float>\n"
        "rho = <float>\n"
        "uinf = <float>\n"
        "re = <float>\n"
        "ret = <float>\n"
        "dt = <float>\n"
        "g_limit0 = <float>\n"
        "g_limit1 = <float>\n"
        "split_t = <integer>\n"
        "zplus = <list of integers>\n"
        "chplot = <string>\n"
        "in_path = <string>\n"
        "out_path = <string>\n"
        "rans_path = <string>\n\n"
        "Example:\n"
        "xlen = 1.0\n"
        "ylen = 2.0\n"
        "zlen = 3.0\n"
        "rho = 1.225\n"
        "uinf = 0.5\n"
        "re = 10000\n"
        "ret = 200\n"
        "dt = 0.01\n"
        "glimit0 = 10.3\n"
        "glimit1 = 20.5\n"
        "split_t = 10\n"
        "zplus = [10, 20, 30]\n"
        "chplot = 'normal' or 'all'\n"
        "in_path = /path/to/input_directory\n"
        "out_path = /path/to/output_directory\n"
        "rans_path = /path/to/rans_datas_file"
    )
    messagebox.showinfo("Parameter File Format", format_text)
    
def show_parameters():
    # Display the contents of parameters.dat file
    
    base_path = os.path.join(os.path.dirname(__file__), "../parameters/")
    
    try:
        with open(os.path.join(base_path, "parameters.dat"), "r") as f:
            params = f.read()
        messagebox.showinfo("Current Parameters", params)
    except FileNotFoundError:
        messagebox.showerror("Error", "parameters.dat file not found!")
    except Exception as e:
        messagebox.showerror("Error", str(e))
        
    
        

def apply_model():
    selected_model = model_var.get()

    # Base path to the directory containing the C files
    base_path = os.path.join(os.path.dirname(__file__), "../parameters/")

    if selected_model == "WRLES_Retau395":
        source_file = os.path.join(base_path, "write_parameters_wrles_395.c")
        output_file = os.path.join(base_path, "write_parameters_wrles_395")
    elif selected_model == "WRLES_Retau1000":
        source_file = os.path.join(base_path, "write_parameters_wrles_1000.c")
        output_file = os.path.join(base_path, "write_parameters_wrles_1000")
    elif selected_model == "WMLES_Retau1000":
        source_file = os.path.join(base_path, "write_parameters_wmles_1000.c")
        output_file = os.path.join(base_path, "write_parameters_wmles_1000")
    else:
        messagebox.showerror("Error", "Please select a model.")
        return

    try:
        # Compile the C file
        compile_command = ["gcc", "-o", f"{output_file}", f"{source_file}"]
        result = subprocess.run(compile_command, capture_output=True, text=True)
        
        if result.returncode != 0:
            error_message = result.stderr
            print(f"Compilation error:\n{error_message}\n")
            messagebox.showerror("Compilation Error", error_message)
            return
        
        print(f"Compiled {source_file} successfully.\n")

        # Run the compiled executable
        result = subprocess.run([output_file], cwd=base_path, capture_output=True, text=True)
        
        
        if result.returncode != 0:
            error_message = result.stderr
            print(f"Runtime error:\n{error_message}\n")
            messagebox.showerror("Runtime Error", error_message)
            return
        
        with open(os.path.join(base_path, "parameters.dat"), "a") as file:
            file.write(f"split_time = {split_time_var.get()}\n")
            file.write(f"nb_split_t = {n_var.get()}\n")
            file.write(f'chplot = {chplot_var.get()}\n')
            file.close()
        
        messagebox.showinfo("Success", f"{selected_model} parameters applied successfully!")
        print(f"Ran {output_file} successfully.\n")
        
    except Exception as e:
        messagebox.showerror("Error", str(e))
        print(f"Error: {str(e)}\n")
        


# Create the main window
root = tk.Tk()
root.title("Run Simulation")

# Styling variables
font = ("new century schoolbook", 13)
bg_color = "#e6f2ff"  # Light blue background
button_color = "#4da6ff"  # Blue button
text_color = "#003366"  # Dark blue text

# Set the window background color
root.configure(bg=bg_color)

# Create and place labels and buttons
tk.Label(root, text="Cflow analysis with the 'parameters.dat' file", font=font, bg=bg_color, fg=text_color).grid(row=0, column=0, columnspan = 2, padx=10, pady=10)

# Dropdown menu for model selection
model_var = tk.StringVar(root)
model_var.set("Select Model")  # Default value

model_menu = tk.OptionMenu(root, model_var, "WRLES_Retau395", "WRLES_Retau1000", "WMLES_Retau1000")
model_menu.config(font=font)
model_menu.grid(row=1, column=0, padx=5, pady=10)

# Frame of split_time parameters
frame_split_time = tk.Frame(root, bg=bg_color)
frame_split_time.grid(row=1, column=1, padx=10, pady=10)

time_split_info = tk.Label(frame_split_time, text="Do you want to split the time series in order to use\nless RAM and to reduce the computation time.\nIf 'Yes' choose the number of time step in the splitted\ntime series (default: Yes).", font=font, bg=bg_color)
time_split_info.grid(row=0, column=0, columnspan=2, padx=5, pady=5)

split_time_var = tk.StringVar(root)
split_time_var.set("Y")  # Default value

yes_radiobutton = tk.Radiobutton(frame_split_time, text="Yes", variable=split_time_var, value="Y", font=font)
yes_radiobutton.grid(row=1, column=0, padx=5, pady=5)
no_radiobutton = tk.Radiobutton(frame_split_time, text="No", variable=split_time_var, value="n", font=font)
no_radiobutton.grid(row=2, column=0, padx=5, pady=5)


n_var = tk.IntVar(root)
n_var.set("Select number")
split_time_menu = tk.OptionMenu(frame_split_time, n_var, int(2**10), int(2**11), int(2**12), int(2**13), int(2**14))
split_time_menu.config(font=font)
split_time_menu.grid(row=1, column=1, rowspan=2, padx=5, pady=5)


# chplot choice
chplot_frame = tk.Frame(root, bg=bg_color)
chplot_frame.grid(row=2, column=0, padx=5, pady=5)

text_chplot = tk.Label(chplot_frame, text="Choose if you want to make the computation\nfor allplanes or only 4 plans equally\ndistributed (default: normal): ", bg=bg_color, font=font)
text_chplot.grid(row=0, column=0, columnspan=2, padx=5, pady=5)

chplot_var = tk.StringVar(root)
chplot_var.set("normal")  # Default value

all_radiobutton = tk.Radiobutton(chplot_frame, text="all", variable=chplot_var, value="all", font=font)
all_radiobutton.grid(row=1, column=0, padx=5, pady=5)
normal_radiobutton = tk.Radiobutton(chplot_frame, text="normal", variable=chplot_var, value="normal", font=font)
normal_radiobutton.grid(row=1, column=1, padx=5, pady=5)


# Apply button
button_apply = tk.Button(root, text="Apply", font=font, bg=button_color, fg="white", command=apply_model)
button_apply.grid(row=2, column=1, pady=5)

separator = ttk.Separator(root, orient='horizontal')
separator.grid(row=3, column=0, columnspan=2, padx=10, pady=10, sticky="ew")


# Show file format button
button_format = tk.Button(root, text="Show Data File Format", command=show_format, font=font)
button_format.grid(row=4, column=0, padx=5, pady=10)

# Show parameters button
button_show_params = tk.Button(root, text="Show Parameters File", command=show_parameters, font=font)
button_show_params.grid(row=4, column=1, padx=5, pady=10)

separator = ttk.Separator(root, orient='horizontal')
separator.grid(row=5, column=0, columnspan=2, padx=10, pady=10, sticky="ew")

# Run button
button_run = tk.Button(root, text="Run", font=font, bg=button_color, fg="white", command=start_script)
button_run.grid(row=6, column=0, pady=10, padx=5)


# Run the application
root.mainloop()
