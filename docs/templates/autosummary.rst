{{ fullname | extract_module_docstring | escape | underline}}

.. currentmodule:: {{ module }}

.. automodule:: {{ fullname }}

  {% block modules %}
  {% if modules %}
  .. autosummary::
    :toctree:
    :template: autosummary.rst
  {% for item in modules %}
    {{ item }}
  {%- endfor %}
  {% endif %}
  {% endblock %}