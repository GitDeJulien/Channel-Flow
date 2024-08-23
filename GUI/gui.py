import tkinter as tk
from tkinter import font
from tkinter import messagebox, scrolledtext
from subprocess import Popen, PIPE, STDOUT
import subprocess
import threading
import os
import signal
import pexpect

def run_script():
    # Check if the parameters.dat file exists and execute the shell script
    try:
            
        display_terminal_output("Running the script...\n")

        # Run the shell script in a separate thread to keep the UI responsive
        thread = threading.Thread(target=execute_script)
        thread.start()
        
    except FileNotFoundError:
        messagebox.showerror("Error", "parameters.dat file not found!")
    except Exception as e:
        messagebox.showerror("Error", str(e))
        
def execute_script():
    # Function to execute the shell script and capture its output
    
    base_path = os.path.join(os.path.dirname(__file__), "../stuff/")
    
    process = subprocess.Popen(
        ["bash", os.path.join(base_path, "run.sh")],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    def read_output(stream):
        for line in iter(stream.readline, ''):
            display_terminal_output(line)
        stream.close()

    # Create threads to read stdout and stderr
    stdout_thread = threading.Thread(target=read_output, args=(process.stdout,))
    stderr_thread = threading.Thread(target=read_output, args=(process.stderr,))

    stdout_thread.start()
    stderr_thread.start()

    stdout_thread.join()
    stderr_thread.join()

    process.wait()

    display_terminal_output("Script finished.\n")

def display_terminal_output(message):
    terminal_output.insert(tk.END, message)
    terminal_output.see(tk.END)  # Scroll to the end of the text box
        

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
        "tStart = <float>\n"
        "dt = <float>\n"
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
        "tStart = 0.0\n"
        "dt = 0.01\n"
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
            display_terminal_output(f"Compilation error:\n{error_message}\n")
            messagebox.showerror("Compilation Error", error_message)
            return
        
        display_terminal_output(f"Compiled {source_file} successfully.\n")

        # Run the compiled executable
        result = subprocess.run([output_file], cwd=base_path, capture_output=True, text=True)
        
        if result.returncode != 0:
            error_message = result.stderr
            display_terminal_output(f"Runtime error:\n{error_message}\n")
            messagebox.showerror("Runtime Error", error_message)
            return
        
        messagebox.showinfo("Success", f"{selected_model} parameters applied successfully!")
        display_terminal_output(f"Ran {output_file} successfully.\n")
        
    except Exception as e:
        messagebox.showerror("Error", str(e))
        display_terminal_output(f"Error: {str(e)}\n")



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
tk.Label(root, text="Cflow analysis with the 'parameters.dat' file", font=font, bg=bg_color, fg=text_color).pack(pady=20)

# Dropdown menu for model selection
model_var = tk.StringVar(root)
model_var.set("Select Model")  # Default value

model_menu = tk.OptionMenu(root, model_var, "WRLES_Retau395", "WRLES_Retau1000", "WMLES_Retau1000")
model_menu.config(font=font)
model_menu.pack(pady=10)

# Apply button
button_apply = tk.Button(root, text="Apply", font=font, bg=button_color, fg="white", command=apply_model)
button_apply.pack(pady=10)

# Show file format button
button_format = tk.Button(root, text="Show Data File Format", command=show_format, font=font)
button_format.pack(pady=10)

# Show parameters button
button_show_params = tk.Button(root, text="Show Parameters File", command=show_parameters, font=font)
button_show_params.pack(pady=10)

# Run button
button_run = tk.Button(root, text="Run", font=font, bg=button_color, fg="white", command=run_script)
button_run.pack(pady=10)

# Add a ScrolledText widget to display terminal output
terminal_output = scrolledtext.ScrolledText(root, height=15, font=font, bg="#1e1e1e", fg="#c5c5c5", insertbackground="white")
terminal_output.pack(pady=20, padx=10, fill=tk.BOTH, expand=True)



# Run the application
root.mainloop()
