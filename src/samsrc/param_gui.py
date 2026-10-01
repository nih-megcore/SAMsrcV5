"""Tkinter editor for SAM analysis parameter files."""

from __future__ import annotations

import os
from pathlib import Path
import sys
import tempfile

from .param_editor import (
    CATEGORIES,
    PROGRAMS,
    SPECS,
    ParameterDocument,
    changed_updates,
    ensure_save_root,
    managed_parameter_path,
    serialize_new,
    validate_values,
)


def _run_gui(tk, ttk, filedialog, messagebox) -> None:
    class ParameterGUI:
        def __init__(self, root, save_root: Path) -> None:
            self.root = root
            self.save_root = save_root
            self.document = ParameterDocument([])
            self.original_values: dict[str, list[str]] = {}
            self.source_path: Path | None = None
            self.fields: dict[str, dict[str, object]] = {}
            self.preview_after: str | None = None

            root.title("SAM Parameter Editor")
            root.geometry("1050x780")
            root.minsize(800, 600)
            self._build_menu()
            self._build_header()
            self._build_editor()
            self._new()

        def _build_menu(self) -> None:
            menu = tk.Menu(self.root)
            file_menu = tk.Menu(menu, tearoff=False)
            file_menu.add_command(label="New", command=self._new, accelerator="Ctrl+N")
            file_menu.add_command(label="Open…", command=self._open, accelerator="Ctrl+O")
            file_menu.add_separator()
            file_menu.add_command(label="Save", command=self._save, accelerator="Ctrl+S")
            file_menu.add_command(label="Save As…", command=self._save_as)
            file_menu.add_separator()
            file_menu.add_command(label="Quit", command=self.root.destroy)
            menu.add_cascade(label="File", menu=file_menu)
            help_menu = tk.Menu(menu, tearoff=False)
            help_menu.add_command(label="About", command=self._about)
            menu.add_cascade(label="Help", menu=help_menu)
            self.root.configure(menu=menu)
            self.root.bind("<Control-n>", lambda _event: self._new())
            self.root.bind("<Control-o>", lambda _event: self._open())
            self.root.bind("<Control-s>", lambda _event: self._save())

        def _build_header(self) -> None:
            header = ttk.Frame(self.root, padding=(10, 8))
            header.pack(fill="x")
            ttk.Label(header, text="Profile:").grid(row=0, column=0, sticky="w")
            self.profile = tk.StringVar(value="sam_cov")
            profile_box = ttk.Combobox(
                header, textvariable=self.profile, values=PROGRAMS, state="readonly", width=12
            )
            profile_box.grid(row=0, column=1, padx=(5, 18), sticky="w")
            profile_box.bind("<<ComboboxSelected>>", self._profile_changed)
            ttk.Label(header, text="Filename:").grid(row=0, column=2, sticky="w")
            self.filename = tk.StringVar(value="analysis.param")
            filename_entry = ttk.Entry(header, textvariable=self.filename, width=34)
            filename_entry.grid(row=0, column=3, padx=5, sticky="ew")
            self.filename.trace_add("write", lambda *_args: self._schedule_preview())
            ttk.Button(header, text="Save", command=self._save).grid(row=0, column=4, padx=(8, 0))
            header.columnconfigure(3, weight=1)
            self.location = ttk.Label(header, text=str(self.save_root), foreground="#555555")
            self.location.grid(row=1, column=0, columnspan=5, pady=(5, 0), sticky="w")

        def _build_editor(self) -> None:
            pane = ttk.Panedwindow(self.root, orient="vertical")
            pane.pack(fill="both", expand=True, padx=10, pady=(0, 10))
            self.notebook = ttk.Notebook(pane)
            pane.add(self.notebook, weight=4)
            category_frames = {}
            for category in CATEGORIES:
                outer = ttk.Frame(self.notebook)
                canvas = tk.Canvas(outer, highlightthickness=0)
                scrollbar = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
                inner = ttk.Frame(canvas, padding=8)
                window = canvas.create_window((0, 0), window=inner, anchor="nw")
                inner.bind(
                    "<Configure>",
                    lambda _event, current=canvas: current.configure(scrollregion=current.bbox("all")),
                )
                canvas.bind(
                    "<Configure>",
                    lambda event, current=canvas, item=window: current.itemconfigure(item, width=event.width),
                )
                canvas.configure(yscrollcommand=scrollbar.set)
                canvas.pack(side="left", fill="both", expand=True)
                scrollbar.pack(side="right", fill="y")
                self.notebook.add(outer, text=category)
                category_frames[category] = inner

            row_by_category = {category: 0 for category in CATEGORIES}
            for spec in SPECS:
                parent = category_frames[spec.category]
                row = row_by_category[spec.category]
                row_by_category[spec.category] += 2
                self._add_field(parent, row, spec)

            preview_frame = ttk.Frame(pane, padding=(0, 8, 0, 0))
            pane.add(preview_frame, weight=2)
            ttk.Label(preview_frame, text="Parameter file preview").pack(anchor="w")
            self.preview = tk.Text(preview_frame, height=12, wrap="none", state="disabled")
            preview_scroll = ttk.Scrollbar(preview_frame, orient="vertical", command=self.preview.yview)
            self.preview.configure(yscrollcommand=preview_scroll.set)
            self.preview.pack(side="left", fill="both", expand=True, pady=(4, 0))
            preview_scroll.pack(side="right", fill="y", pady=(4, 0))
            self.status = ttk.Label(self.root, padding=(10, 0, 10, 8))
            self.status.pack(fill="x")

        def _add_field(self, parent, row: int, spec) -> None:
            enabled = tk.BooleanVar(value=False)
            state: dict[str, object] = {"enabled": enabled, "spec": spec}
            if spec.kind == "flag":
                control = ttk.Checkbutton(parent, text=spec.key, variable=enabled)
                control.grid(row=row, column=0, columnspan=2, sticky="w", pady=4)
                state["label"] = control
            else:
                include = ttk.Checkbutton(parent, variable=enabled)
                include.grid(row=row, column=0, sticky="nw", pady=5)
                label = ttk.Label(parent, text=spec.key, width=22)
                label.grid(row=row, column=1, sticky="nw", pady=6)
                state["label"] = label
                if spec.kind == "repeat":
                    value = tk.Text(parent, height=3, width=55, wrap="none")
                    value.grid(row=row, column=2, sticky="ew", pady=3)
                    value.bind("<<Modified>>", lambda event: self._text_changed(event.widget))
                    state["text"] = value
                else:
                    variable = tk.StringVar()
                    state["value"] = variable
                    if spec.kind == "choice":
                        value = ttk.Combobox(parent, textvariable=variable, values=spec.choices, state="readonly")
                    elif spec.kind == "model":
                        value = ttk.Combobox(
                            parent,
                            textvariable=variable,
                            values=("Nolte 16", "MultiSphere", "SingleSphere 0 0 4"),
                        )
                    elif spec.kind == "imageformat":
                        value = ttk.Combobox(parent, textvariable=variable, values=("ORIG", "TLRC 2"))
                    elif spec.kind == "mu":
                        value = ttk.Combobox(parent, textvariable=variable, values=("+5", "*1"))
                    elif spec.kind == "imagemetric":
                        value = ttk.Combobox(parent, textvariable=variable, values=("Power", "Signal"), state="readonly")
                    else:
                        value = ttk.Entry(parent, textvariable=variable)
                    value.grid(row=row, column=2, sticky="ew", pady=3)
                    if spec.browse:
                        ttk.Button(
                            parent,
                            text="Browse…",
                            command=lambda current=spec, target=variable: self._browse(current, target),
                        ).grid(row=row, column=3, padx=(6, 0), pady=3)
                    variable.trace_add("write", lambda *_args: self._schedule_preview())
                hint = spec.help
                if spec.placeholder:
                    hint = f"{hint}  Format: {spec.placeholder}."
                ttk.Label(parent, text=hint, foreground="#555555", wraplength=650).grid(
                    row=row + 1, column=2, columnspan=2, sticky="w", pady=(0, 5)
                )
            enabled.trace_add("write", lambda *_args: self._schedule_preview())
            self.fields[spec.key] = state
            parent.columnconfigure(2, weight=1)

        def _text_changed(self, widget) -> None:
            if widget.edit_modified():
                widget.edit_modified(False)
                self._schedule_preview()

        def _browse(self, spec, variable) -> None:
            if spec.browse == "directory":
                selected = filedialog.askdirectory(initialdir=variable.get() or str(Path.home()))
            else:
                selected = filedialog.askopenfilename(initialdir=str(Path.home()))
            if selected:
                variable.set(selected)
                self.fields[spec.key]["enabled"].set(True)

        def _profile_changed(self, _event=None) -> None:
            profile = self.profile.get()
            required = {
                "sam_cov": {"CovBand"},
                "sam_wts": {"CovBand", "Model"},
                "sam_3d": {"CovType", "CovBand", "ImageBand", "ImageMetric"},
                "sam_ers": {"Marker", "CovType", "CovBand", "ImageBand", "SmoothBand", "TimeStep", "ImageMetric"},
            }[profile]
            for state in self.fields.values():
                spec = state["spec"]
                label = state["label"]
                if spec.key in required:
                    label.configure(style="Required.TLabel" if spec.kind != "flag" else "Required.TCheckbutton")
                elif profile in spec.programs:
                    label.configure(style="TLabel" if spec.kind != "flag" else "TCheckbutton")
                else:
                    label.configure(style="Irrelevant.TLabel" if spec.kind != "flag" else "Irrelevant.TCheckbutton")
            metric_value = self.fields["ImageMetric"]["value"].get().strip()
            if profile == "sam_3d" and not metric_value:
                state = self.fields["ImageMetric"]
                state["value"].set("Power")
            self._schedule_preview()

        def _values(self) -> dict[str, list[str]]:
            values: dict[str, list[str]] = {}
            for key, state in self.fields.items():
                if not state["enabled"].get():
                    continue
                if "text" in state:
                    raw = state["text"].get("1.0", "end-1c")
                    entries = [line.strip() for line in raw.splitlines() if line.strip()]
                elif state["spec"].kind == "flag":
                    entries = [""]
                else:
                    entries = [state["value"].get().strip()]
                values[key] = entries
            return values

        def _render(self) -> tuple[str, list[str], list[str]]:
            values = self._values()
            errors, warnings = validate_values(values, self.profile.get())
            if self.document.lines:
                updates = changed_updates(self.original_values, values)
                text = self.document.render(updates)
            else:
                text = serialize_new(values)
            return text, errors, warnings

        def _schedule_preview(self) -> None:
            if self.preview_after is not None:
                self.root.after_cancel(self.preview_after)
            self.preview_after = self.root.after(100, self._update_preview)

        def _update_preview(self) -> None:
            self.preview_after = None
            text, errors, warnings = self._render()
            self.preview.configure(state="normal")
            self.preview.delete("1.0", "end")
            self.preview.insert("1.0", text)
            self.preview.configure(state="disabled")
            if errors:
                self.status.configure(text=f"{len(errors)} error(s): {errors[0]}", foreground="#a00000")
            elif warnings:
                self.status.configure(text=f"{len(warnings)} profile warning(s): {warnings[0]}", foreground="#8a5a00")
            else:
                self.status.configure(text="Parameter values are valid for the selected profile.", foreground="#207020")

        def _clear_fields(self) -> None:
            for state in self.fields.values():
                state["enabled"].set(False)
                if "text" in state:
                    state["text"].delete("1.0", "end")
                    state["text"].edit_modified(False)
                if "value" in state:
                    state["value"].set("")

        def _new(self) -> None:
            self.document = ParameterDocument([])
            self.original_values = {}
            self.source_path = None
            self._clear_fields()
            self.filename.set("analysis.param")
            self.profile.set("sam_cov")
            self._profile_changed()

        def _open(self) -> None:
            selected = filedialog.askopenfilename(
                title="Open SAM parameter file",
                initialdir=str(self.save_root),
                filetypes=(("SAM parameter files", "*.param"), ("All files", "*")),
            )
            if not selected:
                return
            path = Path(selected)
            try:
                text = path.read_text(encoding="utf-8")
            except OSError as error:
                messagebox.showerror("Open failed", str(error), parent=self.root)
                return
            self.document = ParameterDocument.parse(text)
            self.original_values = self.document.values()
            self.source_path = path
            self._clear_fields()
            for key, entries in self.original_values.items():
                state = self.fields.get(key)
                if state is None:
                    continue
                state["enabled"].set(True)
                if "text" in state:
                    state["text"].insert("1.0", "\n".join(entries))
                    state["text"].edit_modified(False)
                elif state["spec"].kind != "flag":
                    state["value"].set(entries[0] if entries else "")
            self.filename.set(path.name if path.suffix == ".param" else f"{path.name}.param")
            self.location.configure(text=f"Opened {path}; saves go to {self.save_root}")
            self._profile_changed()

        def _save(self) -> None:
            try:
                path = managed_parameter_path(self.save_root, self.filename.get())
            except ValueError as error:
                messagebox.showerror("Invalid filename", str(error), parent=self.root)
                return
            self._write(path, confirm=path.exists())

        def _save_as(self) -> None:
            selected = filedialog.asksaveasfilename(
                title="Save SAM parameter file",
                initialdir=str(self.save_root),
                initialfile=managed_parameter_path(self.save_root, self.filename.get()).name,
                defaultextension=".param",
                filetypes=(("SAM parameter files", "*.param"),),
            )
            if not selected:
                return
            path = managed_parameter_path(self.save_root, Path(selected).name)
            self.filename.set(path.name)
            self._write(path, confirm=path.exists())

        def _write(self, path: Path, *, confirm: bool) -> None:
            text, errors, warnings = self._render()
            if errors:
                messagebox.showerror("Invalid parameters", "\n".join(errors), parent=self.root)
                return
            if warnings and not messagebox.askyesno(
                "Incomplete profile",
                "The selected profile has warnings:\n\n"
                + "\n".join(warnings)
                + "\n\nSave this partial parameter file?",
                parent=self.root,
            ):
                return
            if confirm and not messagebox.askyesno(
                "Replace parameter file", f"Replace {path}?", parent=self.root
            ):
                return
            temporary: Path | None = None
            try:
                with tempfile.NamedTemporaryFile(
                    mode="w", encoding="utf-8", dir=self.save_root, prefix=".samparam-", delete=False
                ) as stream:
                    stream.write(text)
                    temporary = Path(stream.name)
                os.replace(temporary, path)
            except OSError as error:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
                messagebox.showerror("Save failed", str(error), parent=self.root)
                return
            self.document = ParameterDocument.parse(text)
            self.original_values = self.document.values()
            self.source_path = path
            self.location.configure(text=f"Saved {path}")
            self.status.configure(text=f"Saved {path}", foreground="#207020")

        def _about(self) -> None:
            messagebox.showinfo(
                "About SAM Parameter Editor",
                "Build parameter files for sam_cov, sam_wts, sam_3d, and sam_ers.\n\n"
                "Files are saved in " + str(self.save_root),
                parent=self.root,
            )

    root = tk.Tk()
    style = ttk.Style(root)
    style.configure("Required.TLabel", foreground="#9a2600", font=("TkDefaultFont", 9, "bold"))
    style.configure("Irrelevant.TLabel", foreground="#888888")
    style.configure("Required.TCheckbutton", foreground="#9a2600")
    style.configure("Irrelevant.TCheckbutton", foreground="#888888")
    try:
        save_root = ensure_save_root()
    except OSError as error:
        save_root = Path.home() / "samparams"
        messagebox.showerror("SAM Parameter Editor", f"Cannot create {save_root}:\n{error}", parent=root)
        root.destroy()
        return
    ParameterGUI(root, save_root)
    root.mainloop()


def main() -> None:
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk
    except ImportError as error:
        print(
            "sam_param_gui requires Tk support in the host Python installation "
            "(for example, the operating system's python3-tk package).",
            file=sys.stderr,
        )
        raise SystemExit(1) from error
    try:
        _run_gui(tk, ttk, filedialog, messagebox)
    except tk.TclError as error:
        print(f"sam_param_gui could not open a display: {error}", file=sys.stderr)
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
