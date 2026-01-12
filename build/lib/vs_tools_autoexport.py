import importlib
import inspect
import pkgutil

def export_functions(package_globals):
    """Export functions defined in package modules into the package namespace."""
    package_name = package_globals.get("__name__")
    package_path = package_globals.get("__path__")
    if not package_name or not package_path:
        return

    exported = set(package_globals.get("__all__", []))

    for module_info in pkgutil.iter_modules(package_path):
        if module_info.ispkg or module_info.name.startswith("_"):
            continue
        module = importlib.import_module(f"{package_name}.{module_info.name}")
        for name, obj in vars(module).items():
            if name.startswith("_"):
                continue
            if not inspect.isfunction(obj):
                continue
            if getattr(obj, "__module__", None) != module.__name__:
                continue
            if name in package_globals:
                continue
            package_globals[name] = obj
            exported.add(name)

    package_globals["__all__"] = sorted(exported)
